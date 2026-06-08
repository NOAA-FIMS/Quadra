#!/usr/bin/env bash
set -euo pipefail

target="core/laplace/hessian_structure.hpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  echo "Run install_hessian_structure_dispatch_v1.sh first."
  exit 1
fi

mkdir -p .quadra_patch_backups
cp "$target" ".quadra_patch_backups/hessian_structure.hpp.true_banded_ldlt.$(date +%Y%m%d_%H%M%S).bak"

python3 - <<'PYEOF'
from pathlib import Path

p = Path("core/laplace/hessian_structure.hpp")
s = p.read_text()

old = '''// For now this uses Eigen dense LDLT on a compact dense copy. This is not a
// true O(n*bw^2) banded Cholesky yet. The dispatcher is still useful because it
// provides the location to plug in a true banded backend.
inline double LogDetBandedDenseLDLT(const Eigen::SparseMatrix<double>& H,
                                    const int /*bandwidth*/) {
  Eigen::MatrixXd dense = Eigen::MatrixXd(H);
  Eigen::LDLT<Eigen::MatrixXd> ldlt(dense);

  if (ldlt.info() != Eigen::Success) {
    throw std::runtime_error("Banded dense LDLT failed");
  }

  const auto& D = ldlt.vectorD();
  double logdet = 0.0;
  for (int i = 0; i < D.size(); ++i) {
    if (!(D[i] > 0.0)) {
      throw std::runtime_error("Banded Hessian is not positive definite");
    }
    logdet += std::log(D[i]);
  }

  return logdet;
}
'''

new = '''// True banded LDLT log determinant for symmetric positive definite matrices.
// Uses compact lower-band storage and costs O(n * bandwidth^2).
inline double LogDetBandedLDLT(const Eigen::SparseMatrix<double>& H,
                               const int bandwidth) {
  if (H.rows() != H.cols()) {
    throw std::invalid_argument("Banded LDLT logdet requires square matrix");
  }
  if (bandwidth < 0) {
    throw std::invalid_argument("Banded LDLT bandwidth must be non-negative");
  }

  const int n = static_cast<int>(H.rows());
  if (n == 0) return 0.0;

  const int bw = std::min(bandwidth, n - 1);
  const int stride = bw + 1;

  // Lower-band values: A_band[i * stride + d] = A(i, i - d).
  std::vector<double> A_band(static_cast<std::size_t>(n * stride), 0.0);

  auto at = [&](const int i, const int j) -> double& {
    return A_band[static_cast<std::size_t>(i * stride + (i - j))];
  };

  auto get = [&](const int i, const int j) -> double {
    if (i < j) return 0.0;
    const int d = i - j;
    if (d < 0 || d > bw) return 0.0;
    return A_band[static_cast<std::size_t>(i * stride + d)];
  };

  Eigen::SparseMatrix<double> canonical = H;
  canonical.makeCompressed();

  for (int outer = 0; outer < canonical.outerSize(); ++outer) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(canonical, outer); it; ++it) {
      const int r = static_cast<int>(it.row());
      const int c = static_cast<int>(it.col());
      const double v = it.value();

      const int i = std::max(r, c);
      const int j = std::min(r, c);
      const int d = i - j;

      if (d > bw) {
        if (std::abs(v) > 0.0) {
          throw std::runtime_error("Matrix has nonzero outside declared band");
        }
        continue;
      }

      at(i, j) = v;
    }
  }

  std::vector<double> D(static_cast<std::size_t>(n), 0.0);

  auto L = [&](const int i, const int j) -> double {
    if (i == j) return 1.0;
    if (i < j) return 0.0;
    const int d = i - j;
    if (d <= 0 || d > bw) return 0.0;
    return get(i, j);
  };

  double logdet = 0.0;

  for (int i = 0; i < n; ++i) {
    double diag = get(i, i);

    const int k0_diag = std::max(0, i - bw);
    for (int k = k0_diag; k < i; ++k) {
      const double Lik = L(i, k);
      diag -= Lik * Lik * D[static_cast<std::size_t>(k)];
    }

    if (!(diag > 0.0) || !std::isfinite(diag)) {
      throw std::runtime_error("Banded Hessian is not positive definite");
    }

    D[static_cast<std::size_t>(i)] = diag;
    logdet += std::log(diag);

    const int jmax = std::min(n - 1, i + bw);
    for (int j = i + 1; j <= jmax; ++j) {
      double lij_num = get(j, i);

      const int k0 = std::max(0, std::max(i - bw, j - bw));
      for (int k = k0; k < i; ++k) {
        lij_num -= L(j, k) * D[static_cast<std::size_t>(k)] * L(i, k);
      }

      at(j, i) = lij_num / diag;
    }
  }

  return logdet;
}

// Backwards-compatible name used by early dispatch code.
inline double LogDetBandedDenseLDLT(const Eigen::SparseMatrix<double>& H,
                                    const int bandwidth) {
  return LogDetBandedLDLT(H, bandwidth);
}
'''

if old not in s:
    raise SystemExit("Could not find old banded dense LDLT block")

s = s.replace(old, new, 1)
p.write_text(s)
PYEOF

cat <<'EOF'

Installed true banded LDLT backend.

Run:
  ./run_hessian_structure_dispatch_test.sh

Then rerun no-plus age flat-band benchmark:
  ./run_quadra_flat_band_vs_tmb_age_structured_no_plus_benchmark.sh 10 25,50,100,250,500,1000 10

EOF
