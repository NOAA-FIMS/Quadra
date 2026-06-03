#!/usr/bin/env bash
set -euo pipefail

# install_structure_detector_registry_v1.sh
#
# Adds a reusable structure detector / backend recommendation layer on top of
# the lower-level hessian_structure.hpp work.
#
# New:
#   core/laplace/structure_detector.hpp
#   tests/test_structure_detector_registry.cpp
#   run_structure_detector_registry_test.sh
#
# Purpose:
#   Move from ad hoc model-specific backend choices toward:
#
#     Huu -> StructureDetector -> BackendRecommendation -> Laplace backend
#
# This is intentionally non-invasive. It does not modify existing evaluators yet.

mkdir -p core/laplace tests build/tests

if [[ ! -f core/laplace/hessian_structure.hpp ]]; then
  echo "ERROR: missing core/laplace/hessian_structure.hpp"
  echo "Run install_hessian_structure_dispatch_v1.sh first."
  exit 1
fi

cat > core/laplace/structure_detector.hpp <<'EOF'
#pragma once

#include "hessian_structure.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

namespace quadra {
namespace laplace {

enum class LaplaceBackendKind {
  Diagonal,
  Tridiagonal,
  Banded,
  SparseLDLT,
  DenseLDLT
};

inline const char* ToString(const LaplaceBackendKind backend) {
  switch (backend) {
    case LaplaceBackendKind::Diagonal:
      return "diagonal";
    case LaplaceBackendKind::Tridiagonal:
      return "tridiagonal";
    case LaplaceBackendKind::Banded:
      return "banded";
    case LaplaceBackendKind::SparseLDLT:
      return "sparse_ldlt";
    case LaplaceBackendKind::DenseLDLT:
      return "dense_ldlt";
  }
  return "unknown";
}

struct BackendRecommendation {
  LaplaceBackendKind backend = LaplaceBackendKind::DenseLDLT;
  HessianStructure structure = HessianStructure::Dense;

  int random_size = 0;
  int nnz = 0;
  int bandwidth = 0;
  int max_row_nnz = 0;

  double fill_ratio = 0.0;
  bool symmetric = false;
  bool pattern_reusable = true;
  bool supports_symbolic_reuse = false;
  bool supports_warm_start = true;

  std::string reason;
};

struct StructureDetectorOptions {
  StructureOptions structure_options;

  // Recommended cutoff for using specialized banded code.
  int banded_width_cutoff = 64;

  // Very small random-effect dimensions are usually cheaper as dense even when
  // technically sparse/banded.
  int dense_size_cutoff = 16;

  // Above this fill ratio, dense is usually preferred.
  double dense_fill_ratio = 0.25;

  bool prefer_dense_for_small_matrices = true;
};

class StructureDetector {
 public:
  explicit StructureDetector(
      StructureDetectorOptions options = StructureDetectorOptions())
      : options_(options) {
    options_.structure_options.max_banded_width =
        options_.banded_width_cutoff;
    options_.structure_options.dense_fill_ratio =
        options_.dense_fill_ratio;
  }

  BackendRecommendation Analyze(const Eigen::SparseMatrix<double>& H) const {
    const StructureInfo info =
        InspectHessianStructure(H, options_.structure_options);

    return Recommend(info);
  }

  BackendRecommendation Analyze(const Eigen::MatrixXd& H) const {
    const StructureInfo info =
        InspectHessianStructure(H, options_.structure_options);

    return Recommend(info);
  }

  BackendRecommendation Recommend(const StructureInfo& info) const {
    BackendRecommendation rec;
    rec.random_size = info.rows;
    rec.nnz = info.nnz;
    rec.bandwidth = info.max_bandwidth;
    rec.max_row_nnz = info.max_row_nnz;
    rec.fill_ratio = info.fill_ratio;
    rec.symmetric = info.numerically_symmetric;
    rec.structure = info.detected;

    if (!info.square) {
      rec.backend = LaplaceBackendKind::DenseLDLT;
      rec.pattern_reusable = false;
      rec.supports_symbolic_reuse = false;
      rec.reason = "non-square matrix cannot use Hessian backend";
      return rec;
    }

    if (!info.numerically_symmetric) {
      rec.backend = LaplaceBackendKind::DenseLDLT;
      rec.pattern_reusable = false;
      rec.supports_symbolic_reuse = false;
      rec.reason = "non-symmetric Hessian; falling back to dense";
      return rec;
    }

    if (options_.prefer_dense_for_small_matrices &&
        info.rows <= options_.dense_size_cutoff) {
      rec.backend = LaplaceBackendKind::DenseLDLT;
      rec.supports_symbolic_reuse = false;
      rec.reason = "small matrix; dense LDLT preferred";
      return rec;
    }

    if (info.max_bandwidth == 0) {
      rec.backend = LaplaceBackendKind::Diagonal;
      rec.structure = HessianStructure::Diagonal;
      rec.supports_symbolic_reuse = false;
      rec.reason = "zero off-diagonal bandwidth";
      return rec;
    }

    if (info.max_bandwidth == 1) {
      rec.backend = LaplaceBackendKind::Tridiagonal;
      rec.structure = HessianStructure::Tridiagonal;
      rec.supports_symbolic_reuse = false;
      rec.reason = "unit bandwidth";
      return rec;
    }

    if (info.max_bandwidth <= options_.banded_width_cutoff &&
        info.fill_ratio < options_.dense_fill_ratio) {
      rec.backend = LaplaceBackendKind::Banded;
      rec.structure = HessianStructure::Banded;
      rec.supports_symbolic_reuse = true;
      rec.reason = "fixed narrow band detected";
      return rec;
    }

    if (info.fill_ratio < options_.dense_fill_ratio) {
      rec.backend = LaplaceBackendKind::SparseLDLT;
      rec.structure = HessianStructure::SparsePattern;
      rec.supports_symbolic_reuse = true;
      rec.reason = "general sparse pattern detected";
      return rec;
    }

    rec.backend = LaplaceBackendKind::DenseLDLT;
    rec.structure = HessianStructure::Dense;
    rec.supports_symbolic_reuse = false;
    rec.reason = "high fill ratio; dense LDLT preferred";
    return rec;
  }

  const StructureDetectorOptions& options() const { return options_; }

 private:
  StructureDetectorOptions options_;
};

inline double LogDetWithRecommendation(
    const Eigen::SparseMatrix<double>& H,
    const BackendRecommendation& rec) {
  switch (rec.backend) {
    case LaplaceBackendKind::Diagonal:
      return LogDetDiagonal(H);
    case LaplaceBackendKind::Tridiagonal:
      return LogDetTridiagonalLDLT(H);
    case LaplaceBackendKind::Banded:
      return LogDetBandedLDLT(H, rec.bandwidth);
    case LaplaceBackendKind::SparseLDLT:
      return LogDetSparseLDLT(H);
    case LaplaceBackendKind::DenseLDLT:
      return LogDetDenseLDLT(Eigen::MatrixXd(H));
  }

  throw std::runtime_error("Unknown backend recommendation");
}

}  // namespace laplace
}  // namespace quadra
EOF

cat > tests/test_structure_detector_registry.cpp <<'EOF'
#include "../core/laplace/structure_detector.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using quadra::laplace::BackendRecommendation;
using quadra::laplace::LaplaceBackendKind;
using quadra::laplace::LogDetDenseLDLT;
using quadra::laplace::LogDetWithRecommendation;
using quadra::laplace::StructureDetector;
using quadra::laplace::StructureDetectorOptions;

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

Eigen::SparseMatrix<double> make_sparse(
    int n,
    const std::vector<Eigen::Triplet<double>>& triplets) {
  Eigen::SparseMatrix<double> H(n, n);
  H.setFromTriplets(triplets.begin(), triplets.end());
  return H;
}

StructureDetector make_detector() {
  StructureDetectorOptions opts;
  opts.dense_size_cutoff = 0;  // test structural dispatch even for small cases
  opts.banded_width_cutoff = 8;
  opts.dense_fill_ratio = 0.40;
  opts.prefer_dense_for_small_matrices = false;
  return StructureDetector(opts);
}

void test_diagonal_recommendation() {
  std::vector<Eigen::Triplet<double>> t;
  t.emplace_back(0, 0, 2.0);
  t.emplace_back(1, 1, 3.0);
  t.emplace_back(2, 2, 4.0);

  auto H = make_sparse(3, t);
  auto det = make_detector();
  const auto rec = det.Analyze(H);

  expect_true(rec.backend == LaplaceBackendKind::Diagonal,
              "diagonal recommendation");
  expect_true(rec.bandwidth == 0, "diagonal bandwidth");
  expect_true(rec.pattern_reusable, "diagonal reusable");
  expect_near(LogDetWithRecommendation(H, rec), std::log(24.0), 1e-12,
              "diagonal recommended logdet");
}

void test_tridiagonal_recommendation() {
  const int n = 6;
  std::vector<Eigen::Triplet<double>> t;
  for (int i = 0; i < n; ++i) t.emplace_back(i, i, 4.0);
  for (int i = 0; i < n - 1; ++i) {
    t.emplace_back(i, i + 1, -1.0);
    t.emplace_back(i + 1, i, -1.0);
  }

  auto H = make_sparse(n, t);
  auto det = make_detector();
  const auto rec = det.Analyze(H);

  expect_true(rec.backend == LaplaceBackendKind::Tridiagonal,
              "tridiagonal recommendation");
  expect_true(rec.bandwidth == 1, "tridiagonal bandwidth");

  const double dense = LogDetDenseLDLT(Eigen::MatrixXd(H));
  expect_near(LogDetWithRecommendation(H, rec), dense, 1e-10,
              "tridiagonal recommended logdet");
}

void test_banded_recommendation() {
  const int n = 20;
  std::vector<Eigen::Triplet<double>> t;

  for (int i = 0; i < n; ++i) {
    t.emplace_back(i, i, 8.0);
    for (int d = 1; d <= 4; ++d) {
      if (i + d < n) {
        const double v = -0.05 / static_cast<double>(d);
        t.emplace_back(i, i + d, v);
        t.emplace_back(i + d, i, v);
      }
    }
  }

  auto H = make_sparse(n, t);
  auto det = make_detector();
  const auto rec = det.Analyze(H);

  expect_true(rec.backend == LaplaceBackendKind::Banded,
              "banded recommendation");
  expect_true(rec.bandwidth == 4, "banded bandwidth");
  expect_true(rec.supports_symbolic_reuse, "banded symbolic reuse");

  const double dense = LogDetDenseLDLT(Eigen::MatrixXd(H));
  expect_near(LogDetWithRecommendation(H, rec), dense, 1e-10,
              "banded recommended logdet");
}

void test_sparse_recommendation() {
  const int n = 30;
  std::vector<Eigen::Triplet<double>> t;

  for (int i = 0; i < n; ++i) t.emplace_back(i, i, 5.0);

  // Wide but sparse links.
  t.emplace_back(0, 20, 0.1);
  t.emplace_back(20, 0, 0.1);
  t.emplace_back(3, 27, -0.1);
  t.emplace_back(27, 3, -0.1);

  auto H = make_sparse(n, t);
  auto det = make_detector();
  const auto rec = det.Analyze(H);

  expect_true(rec.backend == LaplaceBackendKind::SparseLDLT,
              "sparse recommendation");
  expect_true(rec.bandwidth > 8, "sparse wide bandwidth");
  expect_true(rec.supports_symbolic_reuse, "sparse symbolic reuse");

  const double dense = LogDetDenseLDLT(Eigen::MatrixXd(H));
  expect_near(LogDetWithRecommendation(H, rec), dense, 1e-10,
              "sparse recommended logdet");
}

void test_dense_recommendation() {
  const int n = 8;
  Eigen::MatrixXd H = Eigen::MatrixXd::Constant(n, n, 0.05);
  H.diagonal().array() = 4.0;

  auto det = make_detector();
  const auto rec = det.Analyze(H);

  expect_true(rec.backend == LaplaceBackendKind::DenseLDLT,
              "dense recommendation");

  std::vector<Eigen::Triplet<double>> t;
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      t.emplace_back(i, j, H(i, j));
    }
  }
  auto S = make_sparse(n, t);

  const double dense = LogDetDenseLDLT(H);
  expect_near(LogDetWithRecommendation(S, rec), dense, 1e-10,
              "dense recommended logdet");
}

void test_small_matrix_dense_preference() {
  StructureDetectorOptions opts;
  opts.prefer_dense_for_small_matrices = true;
  opts.dense_size_cutoff = 16;
  opts.banded_width_cutoff = 8;
  opts.dense_fill_ratio = 0.40;

  StructureDetector det(opts);

  std::vector<Eigen::Triplet<double>> t;
  t.emplace_back(0, 0, 2.0);
  t.emplace_back(1, 1, 3.0);
  t.emplace_back(2, 2, 4.0);

  auto H = make_sparse(3, t);
  const auto rec = det.Analyze(H);

  expect_true(rec.backend == LaplaceBackendKind::DenseLDLT,
              "small matrix dense preference");
}

}  // namespace

int main() {
  test_diagonal_recommendation();
  test_tridiagonal_recommendation();
  test_banded_recommendation();
  test_sparse_recommendation();
  test_dense_recommendation();
  test_small_matrix_dense_preference();

  std::cout << "structure detector registry tests passed\n";
  return 0;
}
EOF

cat > run_structure_detector_registry_test.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/tests

set -x
"${CXX}" ${CXXFLAGS} \
  -Iexternal/Eigen \
  -I. \
  tests/test_structure_detector_registry.cpp \
  -o build/tests/test_structure_detector_registry

./build/tests/test_structure_detector_registry
EOF

chmod +x run_structure_detector_registry_test.sh

cat <<'EOF'

Installed structure detector registry.

Run:
  ./run_structure_detector_registry_test.sh

Next:
  expose BackendRecommendation inside PersistentLaplaceCache adapters,
  then update real Laplace evaluators to call StructureDetector once and reuse the recommendation.

EOF
