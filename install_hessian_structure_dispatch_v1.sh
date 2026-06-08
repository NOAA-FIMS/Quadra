#!/usr/bin/env bash
set -euo pipefail

# install_hessian_structure_dispatch_v1.sh
#
# Adds Quadra's first automatic Hessian structure detector + backend dispatcher.
#
# New file:
#   core/laplace/hessian_structure.hpp
#
# Provides:
#   enum HessianStructure
#   StructureInfo inspect_hessian_structure(...)
#   choose_factorization_backend(...)
#   logdet_diagonal(...)
#   logdet_tridiagonal_ldlt(...)
#   logdet_banded_dense_ldlt(...)
#   logdet_sparse_ldlt(...)
#   logdet_dense_ldlt(...)
#   automatic_logdet(...)
#
# New test:
#   tests/test_hessian_structure_dispatch.cpp
#
# New runner:
#   run_hessian_structure_dispatch_test.sh
#
# This is deliberately standalone and non-invasive. It does not yet modify the
# existing Laplace evaluator. First we prove classification and backend math.

mkdir -p core/laplace tests build/tests

cat > core/laplace/hessian_structure.hpp <<'EOF'
#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>

namespace quadra {
namespace laplace {

enum class HessianStructure {
  Diagonal,
  Tridiagonal,
  Banded,
  SparsePattern,
  Dense
};

inline const char* ToString(const HessianStructure s) {
  switch (s) {
    case HessianStructure::Diagonal:
      return "diagonal";
    case HessianStructure::Tridiagonal:
      return "tridiagonal";
    case HessianStructure::Banded:
      return "banded";
    case HessianStructure::SparsePattern:
      return "sparse_pattern";
    case HessianStructure::Dense:
      return "dense";
  }
  return "unknown";
}

struct StructureInfo {
  int rows = 0;
  int cols = 0;
  int nnz = 0;
  int diagonal_nnz = 0;
  int offdiagonal_nnz = 0;
  int max_bandwidth = 0;
  int max_row_nnz = 0;
  double fill_ratio = 0.0;
  double max_abs_asymmetry = 0.0;
  bool square = false;
  bool structurally_symmetric = true;
  bool numerically_symmetric = true;
  HessianStructure detected = HessianStructure::Dense;
};

struct StructureOptions {
  double zero_tol = 1e-12;
  double symmetry_tol = 1e-10;

  // If max_bandwidth <= this, choose a banded backend.
  int max_banded_width = 64;

  // If fill ratio is above this, dense is usually cheaper/simpler.
  double dense_fill_ratio = 0.25;
};

inline StructureInfo InspectHessianStructure(
    const Eigen::SparseMatrix<double>& H,
    const StructureOptions& options = StructureOptions()) {
  StructureInfo info;
  info.rows = static_cast<int>(H.rows());
  info.cols = static_cast<int>(H.cols());
  info.square = (info.rows == info.cols);

  if (!info.square) {
    throw std::invalid_argument("Hessian structure inspection requires a square matrix");
  }

  Eigen::SparseMatrix<double> canonical = H;
  canonical.makeCompressed();

  int max_row_nnz = 0;
  int current_outer_count = 0;
  int current_outer = -1;

  // For col-major sparse matrices, outer index is column. Row nnz is counted
  // separately below.
  Eigen::VectorXi row_counts = Eigen::VectorXi::Zero(info.rows);

  for (int outer = 0; outer < canonical.outerSize(); ++outer) {
    current_outer_count = 0;

    for (Eigen::SparseMatrix<double>::InnerIterator it(canonical, outer); it; ++it) {
      const int i = static_cast<int>(it.row());
      const int j = static_cast<int>(it.col());
      const double v = it.value();

      if (std::abs(v) <= options.zero_tol) {
        continue;
      }

      ++info.nnz;
      ++current_outer_count;
      ++row_counts[i];

      if (i == j) {
        ++info.diagonal_nnz;
      } else {
        ++info.offdiagonal_nnz;
      }

      info.max_bandwidth = std::max(info.max_bandwidth, std::abs(i - j));
    }

    (void)current_outer;
  }

  for (int i = 0; i < row_counts.size(); ++i) {
    max_row_nnz = std::max(max_row_nnz, row_counts[i]);
  }
  info.max_row_nnz = max_row_nnz;

  const double total = static_cast<double>(info.rows) * static_cast<double>(info.cols);
  info.fill_ratio = total > 0.0 ? static_cast<double>(info.nnz) / total : 0.0;

  // Numeric symmetry check. This is O(nnz * lookup), acceptable for detection
  // and tests. Production evaluators can cache this after pattern discovery.
  info.max_abs_asymmetry = 0.0;
  for (int outer = 0; outer < canonical.outerSize(); ++outer) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(canonical, outer); it; ++it) {
      const int i = static_cast<int>(it.row());
      const int j = static_cast<int>(it.col());
      const double v = it.value();

      if (std::abs(v) <= options.zero_tol) {
        continue;
      }

      const double vt = canonical.coeff(j, i);
      const double diff = std::abs(v - vt);
      info.max_abs_asymmetry = std::max(info.max_abs_asymmetry, diff);

      if (std::abs(vt) <= options.zero_tol) {
        info.structurally_symmetric = false;
      }
    }
  }

  info.numerically_symmetric = info.max_abs_asymmetry <= options.symmetry_tol;

  if (info.max_bandwidth == 0) {
    info.detected = HessianStructure::Diagonal;
  } else if (info.max_bandwidth == 1) {
    info.detected = HessianStructure::Tridiagonal;
  } else if (info.max_bandwidth <= options.max_banded_width &&
             info.fill_ratio < options.dense_fill_ratio) {
    info.detected = HessianStructure::Banded;
  } else if (info.fill_ratio < options.dense_fill_ratio) {
    info.detected = HessianStructure::SparsePattern;
  } else {
    info.detected = HessianStructure::Dense;
  }

  return info;
}

inline StructureInfo InspectHessianStructure(
    const Eigen::MatrixXd& H,
    const StructureOptions& options = StructureOptions()) {
  if (H.rows() != H.cols()) {
    throw std::invalid_argument("Hessian structure inspection requires a square matrix");
  }

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(static_cast<std::size_t>(H.rows() * H.cols()));

  for (int i = 0; i < H.rows(); ++i) {
    for (int j = 0; j < H.cols(); ++j) {
      if (std::abs(H(i, j)) > options.zero_tol) {
        triplets.emplace_back(i, j, H(i, j));
      }
    }
  }

  Eigen::SparseMatrix<double> S(H.rows(), H.cols());
  S.setFromTriplets(triplets.begin(), triplets.end());
  return InspectHessianStructure(S, options);
}

inline HessianStructure ChooseFactorizationBackend(
    const StructureInfo& info,
    const StructureOptions& options = StructureOptions()) {
  if (!info.square) {
    throw std::invalid_argument("Cannot factor non-square Hessian");
  }
  if (!info.numerically_symmetric) {
    throw std::invalid_argument("Cannot use symmetric Hessian backend on non-symmetric matrix");
  }

  if (info.max_bandwidth == 0) {
    return HessianStructure::Diagonal;
  }
  if (info.max_bandwidth == 1) {
    return HessianStructure::Tridiagonal;
  }
  if (info.max_bandwidth <= options.max_banded_width &&
      info.fill_ratio < options.dense_fill_ratio) {
    return HessianStructure::Banded;
  }
  if (info.fill_ratio < options.dense_fill_ratio) {
    return HessianStructure::SparsePattern;
  }
  return HessianStructure::Dense;
}

inline double LogDetDiagonal(const Eigen::SparseMatrix<double>& H) {
  if (H.rows() != H.cols()) {
    throw std::invalid_argument("Diagonal logdet requires square matrix");
  }

  double logdet = 0.0;

  for (int i = 0; i < H.rows(); ++i) {
    const double d = H.coeff(i, i);
    if (!(d > 0.0)) {
      throw std::runtime_error("Diagonal Hessian is not positive definite");
    }
    logdet += std::log(d);
  }

  return logdet;
}

// LDLT logdet for symmetric tridiagonal positive definite matrices.
// This avoids sparse symbolic overhead for max_bandwidth == 1.
inline double LogDetTridiagonalLDLT(const Eigen::SparseMatrix<double>& H) {
  if (H.rows() != H.cols()) {
    throw std::invalid_argument("Tridiagonal logdet requires square matrix");
  }

  const int n = static_cast<int>(H.rows());
  if (n == 0) return 0.0;

  double logdet = 0.0;

  double d_prev = H.coeff(0, 0);
  if (!(d_prev > 0.0)) {
    throw std::runtime_error("Tridiagonal Hessian is not positive definite");
  }
  logdet += std::log(d_prev);

  for (int i = 1; i < n; ++i) {
    const double e = H.coeff(i, i - 1);
    const double diag = H.coeff(i, i);

    const double d = diag - (e * e) / d_prev;
    if (!(d > 0.0)) {
      throw std::runtime_error("Tridiagonal Hessian is not positive definite");
    }

    logdet += std::log(d);
    d_prev = d;
  }

  return logdet;
}

// For now this uses Eigen dense LDLT on a compact dense copy. This is not a
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

inline double LogDetSparseLDLT(const Eigen::SparseMatrix<double>& H) {
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
  ldlt.compute(H);

  if (ldlt.info() != Eigen::Success) {
    throw std::runtime_error("Sparse LDLT failed");
  }

  const auto& D = ldlt.vectorD();
  double logdet = 0.0;
  for (int i = 0; i < D.size(); ++i) {
    if (!(D[i] > 0.0)) {
      throw std::runtime_error("Sparse Hessian is not positive definite");
    }
    logdet += std::log(D[i]);
  }

  return logdet;
}

inline double LogDetDenseLDLT(const Eigen::MatrixXd& H) {
  Eigen::LDLT<Eigen::MatrixXd> ldlt(H);

  if (ldlt.info() != Eigen::Success) {
    throw std::runtime_error("Dense LDLT failed");
  }

  const auto& D = ldlt.vectorD();
  double logdet = 0.0;
  for (int i = 0; i < D.size(); ++i) {
    if (!(D[i] > 0.0)) {
      throw std::runtime_error("Dense Hessian is not positive definite");
    }
    logdet += std::log(D[i]);
  }

  return logdet;
}

inline double AutomaticLogDet(
    const Eigen::SparseMatrix<double>& H,
    const StructureOptions& options = StructureOptions(),
    HessianStructure* selected_backend = nullptr,
    StructureInfo* out_info = nullptr) {
  const StructureInfo info = InspectHessianStructure(H, options);
  const HessianStructure backend = ChooseFactorizationBackend(info, options);

  if (selected_backend != nullptr) {
    *selected_backend = backend;
  }
  if (out_info != nullptr) {
    *out_info = info;
  }

  switch (backend) {
    case HessianStructure::Diagonal:
      return LogDetDiagonal(H);
    case HessianStructure::Tridiagonal:
      return LogDetTridiagonalLDLT(H);
    case HessianStructure::Banded:
      return LogDetBandedDenseLDLT(H, info.max_bandwidth);
    case HessianStructure::SparsePattern:
      return LogDetSparseLDLT(H);
    case HessianStructure::Dense:
      return LogDetDenseLDLT(Eigen::MatrixXd(H));
  }

  throw std::runtime_error("Unknown Hessian backend");
}

}  // namespace laplace
}  // namespace quadra
EOF

cat > tests/test_hessian_structure_dispatch.cpp <<'EOF'
#include "../core/laplace/hessian_structure.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using quadra::laplace::AutomaticLogDet;
using quadra::laplace::ChooseFactorizationBackend;
using quadra::laplace::HessianStructure;
using quadra::laplace::InspectHessianStructure;
using quadra::laplace::StructureInfo;
using quadra::laplace::StructureOptions;
using quadra::laplace::ToString;

namespace {

void expect_true(bool cond, const char* msg) {
  if (!cond) {
    throw std::runtime_error(msg);
  }
}

void expect_near(double a, double b, double tol, const char* msg) {
  if (std::abs(a - b) > tol) {
    std::cerr << msg << ": a=" << a << " b=" << b
              << " diff=" << std::abs(a - b) << "\n";
    throw std::runtime_error(msg);
  }
}

Eigen::SparseMatrix<double> make_sparse_from_triplets(
    int n,
    const std::vector<Eigen::Triplet<double>>& triplets) {
  Eigen::SparseMatrix<double> H(n, n);
  H.setFromTriplets(triplets.begin(), triplets.end());
  return H;
}

void test_diagonal() {
  std::vector<Eigen::Triplet<double>> t;
  t.emplace_back(0, 0, 2.0);
  t.emplace_back(1, 1, 3.0);
  t.emplace_back(2, 2, 4.0);

  auto H = make_sparse_from_triplets(3, t);
  StructureInfo info;
  HessianStructure backend;
  const double logdet = AutomaticLogDet(H, StructureOptions(), &backend, &info);

  expect_true(info.detected == HessianStructure::Diagonal, "diagonal detected");
  expect_true(backend == HessianStructure::Diagonal, "diagonal backend");
  expect_true(info.max_bandwidth == 0, "diagonal bandwidth");
  expect_near(logdet, std::log(24.0), 1e-12, "diagonal logdet");
}

void test_tridiagonal() {
  std::vector<Eigen::Triplet<double>> t;
  const int n = 4;
  for (int i = 0; i < n; ++i) t.emplace_back(i, i, 4.0);
  for (int i = 0; i < n - 1; ++i) {
    t.emplace_back(i, i + 1, -1.0);
    t.emplace_back(i + 1, i, -1.0);
  }

  auto H = make_sparse_from_triplets(n, t);
  StructureInfo info;
  HessianStructure backend;
  const double auto_logdet = AutomaticLogDet(H, StructureOptions(), &backend, &info);

  const double dense_logdet =
      quadra::laplace::LogDetDenseLDLT(Eigen::MatrixXd(H));

  expect_true(info.detected == HessianStructure::Tridiagonal, "tridiagonal detected");
  expect_true(backend == HessianStructure::Tridiagonal, "tridiagonal backend");
  expect_true(info.max_bandwidth == 1, "tridiagonal bandwidth");
  expect_near(auto_logdet, dense_logdet, 1e-10, "tridiagonal logdet");
}

void test_banded() {
  const int n = 10;
  std::vector<Eigen::Triplet<double>> t;

  for (int i = 0; i < n; ++i) {
    t.emplace_back(i, i, 10.0);
    for (int d = 1; d <= 3; ++d) {
      if (i + d < n) {
        const double v = -0.1 / static_cast<double>(d);
        t.emplace_back(i, i + d, v);
        t.emplace_back(i + d, i, v);
      }
    }
  }

  auto H = make_sparse_from_triplets(n, t);
  StructureOptions opts;
  opts.max_banded_width = 8;
  opts.dense_fill_ratio = 0.60;

  StructureInfo info;
  HessianStructure backend;
  const double auto_logdet = AutomaticLogDet(H, opts, &backend, &info);
  const double dense_logdet =
      quadra::laplace::LogDetDenseLDLT(Eigen::MatrixXd(H));

  expect_true(info.detected == HessianStructure::Banded, "banded detected");
  expect_true(backend == HessianStructure::Banded, "banded backend");
  expect_true(info.max_bandwidth == 3, "banded bandwidth");
  expect_near(auto_logdet, dense_logdet, 1e-10, "banded logdet");
}

void test_sparse_pattern() {
  const int n = 20;
  std::vector<Eigen::Triplet<double>> t;
  for (int i = 0; i < n; ++i) t.emplace_back(i, i, 5.0);

  // Wide sparse links, still SPD by diagonal dominance.
  t.emplace_back(0, 10, 0.1);
  t.emplace_back(10, 0, 0.1);
  t.emplace_back(3, 18, -0.1);
  t.emplace_back(18, 3, -0.1);

  auto H = make_sparse_from_triplets(n, t);
  StructureOptions opts;
  opts.max_banded_width = 4;
  opts.dense_fill_ratio = 0.25;

  StructureInfo info;
  HessianStructure backend;
  const double auto_logdet = AutomaticLogDet(H, opts, &backend, &info);
  const double sparse_logdet = quadra::laplace::LogDetSparseLDLT(H);

  expect_true(info.detected == HessianStructure::SparsePattern, "sparse detected");
  expect_true(backend == HessianStructure::SparsePattern, "sparse backend");
  expect_near(auto_logdet, sparse_logdet, 1e-10, "sparse logdet");
}

void test_dense() {
  Eigen::MatrixXd H = Eigen::MatrixXd::Identity(5, 5) * 6.0;
  H.array() += 0.05;
  H = 0.5 * (H + H.transpose());
  H.diagonal().array() += 1.0;

  StructureOptions opts;
  opts.dense_fill_ratio = 0.25;

  const auto info = InspectHessianStructure(H, opts);
  const auto backend = ChooseFactorizationBackend(info, opts);
  const double dense_logdet = quadra::laplace::LogDetDenseLDLT(H);

  std::vector<Eigen::Triplet<double>> t;
  for (int i = 0; i < H.rows(); ++i) {
    for (int j = 0; j < H.cols(); ++j) {
      t.emplace_back(i, j, H(i, j));
    }
  }

  Eigen::SparseMatrix<double> S(H.rows(), H.cols());
  S.setFromTriplets(t.begin(), t.end());

  HessianStructure auto_backend;
  const double auto_logdet = AutomaticLogDet(S, opts, &auto_backend);

  expect_true(info.detected == HessianStructure::Dense, "dense detected");
  expect_true(backend == HessianStructure::Dense, "dense backend");
  expect_true(auto_backend == HessianStructure::Dense, "auto dense backend");
  expect_near(auto_logdet, dense_logdet, 1e-10, "dense logdet");
}

}  // namespace

int main() {
  test_diagonal();
  test_tridiagonal();
  test_banded();
  test_sparse_pattern();
  test_dense();

  std::cout << "hessian structure dispatch tests passed\n";
  return 0;
}
EOF

cat > run_hessian_structure_dispatch_test.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/tests

set -x
"${CXX}" ${CXXFLAGS} \
  -Iexternal/Eigen \
  -I. \
  tests/test_hessian_structure_dispatch.cpp \
  -o build/tests/test_hessian_structure_dispatch

./build/tests/test_hessian_structure_dispatch
EOF

chmod +x run_hessian_structure_dispatch_test.sh

cat <<'EOF'

Installed Hessian structure detection + dispatch scaffold.

Run:
  ./run_hessian_structure_dispatch_test.sh

Next after tests pass:
  wire AutomaticLogDet / backend choice into Laplace evaluators,
  then replace LogDetBandedDenseLDLT with a true banded LDLT backend.

EOF
