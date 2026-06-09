#!/usr/bin/env bash
set -euo pipefail

echo "== Add selected-inverse diagonal utility for random-effect covariance =="

stamp="$(date +%Y%m%d_%H%M%S)"
mkdir -p .quadra_patch_backups core/uncertainty tests build/tests

hdr="core/uncertainty/selected_inverse_diagonal.hpp"
test_src="tests/test_selected_inverse_diagonal.cpp"

[[ -f "$hdr" ]] && cp "$hdr" ".quadra_patch_backups/selected_inverse_diagonal.hpp.${stamp}.bak"
[[ -f "$test_src" ]] && cp "$test_src" ".quadra_patch_backups/test_selected_inverse_diagonal.cpp.${stamp}.bak"

cat > "$hdr" <<'HPP'
#pragma once

#include <Eigen/Core>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>

#include <cmath>
#include <limits>
#include <string>
#include <vector>

namespace quadra {
namespace uncertainty {

struct SelectedInverseDiagonalResult {
  std::vector<double> variance;
  std::vector<double> standard_error;
  bool success = false;
  std::string message;
};

// Level-1 conditional random-effect covariance utility.
//
// For a positive-definite random-effect Hessian H_uu, the conditional
// covariance is approximately inv(H_uu). This conservative implementation
// extracts diag(inv(H_uu)) by solving H_uu x_i = e_i and reading x_i[i].
//
// This is deliberately simple and validation-friendly. Later, this can be
// replaced with a true selected-inverse algorithm for very large models.
inline SelectedInverseDiagonalResult selected_inverse_diagonal_from_spd_hessian(
    const Eigen::SparseMatrix<double>& hessian,
    double min_variance = 0.0) {
  SelectedInverseDiagonalResult out;

  const int n = static_cast<int>(hessian.rows());
  if (hessian.rows() != hessian.cols()) {
    out.message = "Hessian is not square";
    return out;
  }

  out.variance.assign(static_cast<std::size_t>(n),
                      std::numeric_limits<double>::quiet_NaN());
  out.standard_error.assign(static_cast<std::size_t>(n),
                            std::numeric_limits<double>::quiet_NaN());

  if (n == 0) {
    out.success = true;
    out.message = "empty Hessian";
    return out;
  }

  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
  ldlt.analyzePattern(hessian);
  ldlt.factorize(hessian);

  if (ldlt.info() != Eigen::Success) {
    out.message = "SimplicialLDLT factorization failed";
    return out;
  }

  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(n);
  for (int i = 0; i < n; ++i) {
    rhs.setZero();
    rhs[i] = 1.0;

    const Eigen::VectorXd sol = ldlt.solve(rhs);
    if (ldlt.info() != Eigen::Success) {
      out.message = "SimplicialLDLT solve failed";
      return out;
    }

    double v = sol[i];

    if (std::isfinite(v) && v < 0.0 && std::abs(v) <= 1.0e-12) {
      v = 0.0;
    }

    if (std::isfinite(v) && v < min_variance) {
      v = min_variance;
    }

    out.variance[static_cast<std::size_t>(i)] = v;
    out.standard_error[static_cast<std::size_t>(i)] =
        (std::isfinite(v) && v >= 0.0) ? std::sqrt(v)
                                       : std::numeric_limits<double>::quiet_NaN();
  }

  out.success = true;
  out.message = "ok";
  return out;
}

}  // namespace uncertainty
}  // namespace quadra
HPP

cat > "$test_src" <<'CPP'
#include "../core/uncertainty/selected_inverse_diagonal.hpp"

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <cmath>
#include <iostream>
#include <vector>

int main() {
  using Triplet = Eigen::Triplet<double>;

  // SPD tridiagonal matrix:
  // [ 2 -1  0 ]
  // [-1  2 -1 ]
  // [ 0 -1  2 ]
  //
  // inverse diag = [0.75, 1.0, 0.75]
  Eigen::SparseMatrix<double> h(3, 3);
  std::vector<Triplet> triplets;
  triplets.emplace_back(0, 0, 2.0);
  triplets.emplace_back(1, 1, 2.0);
  triplets.emplace_back(2, 2, 2.0);
  triplets.emplace_back(0, 1, -1.0);
  triplets.emplace_back(1, 0, -1.0);
  triplets.emplace_back(1, 2, -1.0);
  triplets.emplace_back(2, 1, -1.0);
  h.setFromTriplets(triplets.begin(), triplets.end());

  const auto result =
      quadra::uncertainty::selected_inverse_diagonal_from_spd_hessian(h);

  if (!result.success) {
    std::cerr << "selected inverse failed: " << result.message << "\n";
    return 1;
  }

  const double expected[] = {0.75, 1.0, 0.75};
  for (int i = 0; i < 3; ++i) {
    const double got = result.variance.at(static_cast<std::size_t>(i));
    if (std::abs(got - expected[i]) > 1.0e-10) {
      std::cerr << "variance mismatch at " << i << ": got " << got
                << ", expected " << expected[i] << "\n";
      return 1;
    }
  }

  std::cout << "selected_inverse_diagonal_test_passed\n";
  std::cout << "variance:";
  for (double v : result.variance) std::cout << " " << v;
  std::cout << "\nstandard_error:";
  for (double se : result.standard_error) std::cout << " " << se;
  std::cout << "\n";
  return 0;
}
CPP

cat > run_selected_inverse_diagonal_test.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/tests

eigen_include=""
for d in external/eigen core/eigen external/eigen3 external/Eigen /opt/homebrew/include/eigen3 /usr/local/include/eigen3 /usr/include/eigen3; do
  if [[ -f "$d/Eigen/Core" ]]; then
    eigen_include="$d"
    break
  fi
done

if [[ -z "$eigen_include" ]]; then
  found="$(find . -path '*/Eigen/Core' -type f 2>/dev/null | head -1 || true)"
  [[ -n "$found" ]] && eigen_include="$(dirname "$(dirname "$found")")"
fi

if [[ -z "$eigen_include" ]]; then
  echo "ERROR: could not find Eigen/Core" >&2
  exit 1
fi

echo "Using Eigen include: $eigen_include"

c++ -std=c++17 -O3 -I"$eigen_include" -I. \
  -o build/tests/test_selected_inverse_diagonal \
  tests/test_selected_inverse_diagonal.cpp

./build/tests/test_selected_inverse_diagonal
SH
chmod +x run_selected_inverse_diagonal_test.sh

cat > inspect_selected_inverse_diagonal_v1.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

echo "== Selected inverse diagonal utility =="
grep -n "selected_inverse_diagonal_from_spd_hessian\\|SelectedInverseDiagonalResult\\|SimplicialLDLT" \
  core/uncertainty/selected_inverse_diagonal.hpp

echo
echo "== Test source =="
grep -n "inverse diag\\|selected_inverse_diagonal_test_passed\\|expected" \
  tests/test_selected_inverse_diagonal.cpp

echo
echo "== Run test =="
./run_selected_inverse_diagonal_test.sh
SH
chmod +x inspect_selected_inverse_diagonal_v1.sh

echo
echo "Installed:"
echo "  $hdr"
echo "  $test_src"
echo "  run_selected_inverse_diagonal_test.sh"
echo "  inspect_selected_inverse_diagonal_v1.sh"
echo
echo "Run:"
echo "  ./inspect_selected_inverse_diagonal_v1.sh"
