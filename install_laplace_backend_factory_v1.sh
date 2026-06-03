#!/usr/bin/env bash
set -euo pipefail

# install_laplace_backend_factory_v1.sh
#
# Adds Quadra's common Laplace backend interface and backend factory.
#
# New:
#   core/laplace/laplace_backend.hpp
#   core/laplace/laplace_backend_factory.hpp
#   tests/test_laplace_backend_factory.cpp
#   run_laplace_backend_factory_test.sh
#
# Depends on:
#   core/laplace/hessian_structure.hpp
#   core/laplace/structure_detector.hpp
#
# This is non-invasive. It does not modify existing evaluators yet.
#
# Goal:
#   One common interface for:
#     diagonal
#     tridiagonal
#     banded
#     sparse LDLT
#     dense LDLT
#
# Then evaluators/cache can store one backend pointer and stop caring which
# structure was selected.

mkdir -p core/laplace tests build/tests

if [[ ! -f core/laplace/hessian_structure.hpp ]]; then
  echo "ERROR: missing core/laplace/hessian_structure.hpp"
  echo "Run install_hessian_structure_dispatch_v1.sh first."
  exit 1
fi

if [[ ! -f core/laplace/structure_detector.hpp ]]; then
  echo "ERROR: missing core/laplace/structure_detector.hpp"
  echo "Run install_structure_detector_registry_v1.sh first."
  exit 1
fi

cat > core/laplace/laplace_backend.hpp <<'EOF'
#pragma once

#include "hessian_structure.hpp"
#include "structure_detector.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace quadra {
namespace laplace {

class LaplaceBackend {
 public:
  virtual ~LaplaceBackend() = default;

  virtual const char* name() const = 0;

  virtual void analyze_pattern(const Eigen::SparseMatrix<double>& H) {
    analyzed_ = true;
    rows_ = static_cast<int>(H.rows());
    cols_ = static_cast<int>(H.cols());
  }

  virtual void factorize(const Eigen::SparseMatrix<double>& H) = 0;

  virtual double logdet() const = 0;

  virtual bool is_spd() const = 0;

  virtual int rows() const { return rows_; }
  virtual int cols() const { return cols_; }
  virtual bool analyzed() const { return analyzed_; }

 protected:
  bool analyzed_ = false;
  int rows_ = 0;
  int cols_ = 0;
};

class DiagonalBackend final : public LaplaceBackend {
 public:
  const char* name() const override { return "diagonal"; }

  void analyze_pattern(const Eigen::SparseMatrix<double>& H) override {
    LaplaceBackend::analyze_pattern(H);
    if (H.rows() != H.cols()) {
      throw std::invalid_argument("DiagonalBackend requires square matrix");
    }
  }

  void factorize(const Eigen::SparseMatrix<double>& H) override {
    analyze_pattern(H);

    logdet_ = 0.0;
    spd_ = true;

    for (int i = 0; i < H.rows(); ++i) {
      const double d = H.coeff(i, i);
      if (!(d > 0.0) || !std::isfinite(d)) {
        spd_ = false;
        logdet_ = std::numeric_limits<double>::quiet_NaN();
        return;
      }
      logdet_ += std::log(d);
    }
  }

  double logdet() const override { return logdet_; }
  bool is_spd() const override { return spd_; }

 private:
  double logdet_ = 0.0;
  bool spd_ = false;
};

class TridiagonalBackend final : public LaplaceBackend {
 public:
  const char* name() const override { return "tridiagonal"; }

  void analyze_pattern(const Eigen::SparseMatrix<double>& H) override {
    LaplaceBackend::analyze_pattern(H);
    if (H.rows() != H.cols()) {
      throw std::invalid_argument("TridiagonalBackend requires square matrix");
    }
  }

  void factorize(const Eigen::SparseMatrix<double>& H) override {
    analyze_pattern(H);

    try {
      logdet_ = LogDetTridiagonalLDLT(H);
      spd_ = true;
    } catch (...) {
      logdet_ = std::numeric_limits<double>::quiet_NaN();
      spd_ = false;
    }
  }

  double logdet() const override { return logdet_; }
  bool is_spd() const override { return spd_; }

 private:
  double logdet_ = 0.0;
  bool spd_ = false;
};

class BandedBackend final : public LaplaceBackend {
 public:
  explicit BandedBackend(int bandwidth) : bandwidth_(bandwidth) {
    if (bandwidth_ < 0) {
      throw std::invalid_argument("BandedBackend bandwidth must be non-negative");
    }
  }

  const char* name() const override { return "banded"; }

  int bandwidth() const { return bandwidth_; }

  void analyze_pattern(const Eigen::SparseMatrix<double>& H) override {
    LaplaceBackend::analyze_pattern(H);
    if (H.rows() != H.cols()) {
      throw std::invalid_argument("BandedBackend requires square matrix");
    }

    // Validate that the pattern fits the declared band.
    for (int outer = 0; outer < H.outerSize(); ++outer) {
      for (Eigen::SparseMatrix<double>::InnerIterator it(H, outer); it; ++it) {
        const int i = static_cast<int>(it.row());
        const int j = static_cast<int>(it.col());
        if (std::abs(i - j) > bandwidth_ && std::abs(it.value()) > 0.0) {
          throw std::invalid_argument("BandedBackend pattern exceeds bandwidth");
        }
      }
    }
  }

  void factorize(const Eigen::SparseMatrix<double>& H) override {
    analyze_pattern(H);

    try {
      logdet_ = LogDetBandedLDLT(H, bandwidth_);
      spd_ = true;
    } catch (...) {
      logdet_ = std::numeric_limits<double>::quiet_NaN();
      spd_ = false;
    }
  }

  double logdet() const override { return logdet_; }
  bool is_spd() const override { return spd_; }

 private:
  int bandwidth_ = 0;
  double logdet_ = 0.0;
  bool spd_ = false;
};

class SparseLDLTBackend final : public LaplaceBackend {
 public:
  const char* name() const override { return "sparse_ldlt"; }

  void analyze_pattern(const Eigen::SparseMatrix<double>& H) override {
    LaplaceBackend::analyze_pattern(H);
    if (H.rows() != H.cols()) {
      throw std::invalid_argument("SparseLDLTBackend requires square matrix");
    }

    Eigen::SparseMatrix<double> canonical = H;
    canonical.makeCompressed();

    ldlt_.analyzePattern(canonical);
    analyzed_ = true;
    symbolic_ready_ = true;
  }

  void factorize(const Eigen::SparseMatrix<double>& H) override {
    Eigen::SparseMatrix<double> canonical = H;
    canonical.makeCompressed();

    if (!symbolic_ready_) {
      analyze_pattern(canonical);
    }

    ldlt_.factorize(canonical);

    if (ldlt_.info() != Eigen::Success) {
      spd_ = false;
      logdet_ = std::numeric_limits<double>::quiet_NaN();
      return;
    }

    const auto& D = ldlt_.vectorD();
    logdet_ = 0.0;
    spd_ = true;

    for (int i = 0; i < D.size(); ++i) {
      if (!(D[i] > 0.0) || !std::isfinite(D[i])) {
        spd_ = false;
        logdet_ = std::numeric_limits<double>::quiet_NaN();
        return;
      }
      logdet_ += std::log(D[i]);
    }
  }

  double logdet() const override { return logdet_; }
  bool is_spd() const override { return spd_; }

 private:
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt_;
  bool symbolic_ready_ = false;
  double logdet_ = 0.0;
  bool spd_ = false;
};

class DenseLDLTBackend final : public LaplaceBackend {
 public:
  const char* name() const override { return "dense_ldlt"; }

  void analyze_pattern(const Eigen::SparseMatrix<double>& H) override {
    LaplaceBackend::analyze_pattern(H);
    if (H.rows() != H.cols()) {
      throw std::invalid_argument("DenseLDLTBackend requires square matrix");
    }
  }

  void factorize(const Eigen::SparseMatrix<double>& H) override {
    analyze_pattern(H);

    Eigen::MatrixXd dense = Eigen::MatrixXd(H);
    Eigen::LDLT<Eigen::MatrixXd> ldlt(dense);

    if (ldlt.info() != Eigen::Success) {
      spd_ = false;
      logdet_ = std::numeric_limits<double>::quiet_NaN();
      return;
    }

    const auto& D = ldlt.vectorD();
    logdet_ = 0.0;
    spd_ = true;

    for (int i = 0; i < D.size(); ++i) {
      if (!(D[i] > 0.0) || !std::isfinite(D[i])) {
        spd_ = false;
        logdet_ = std::numeric_limits<double>::quiet_NaN();
        return;
      }
      logdet_ += std::log(D[i]);
    }
  }

  double logdet() const override { return logdet_; }
  bool is_spd() const override { return spd_; }

 private:
  double logdet_ = 0.0;
  bool spd_ = false;
};

}  // namespace laplace
}  // namespace quadra
EOF

cat > core/laplace/laplace_backend_factory.hpp <<'EOF'
#pragma once

#include "laplace_backend.hpp"
#include "structure_detector.hpp"

#include <memory>
#include <stdexcept>

namespace quadra {
namespace laplace {

inline std::unique_ptr<LaplaceBackend> CreateLaplaceBackend(
    const BackendRecommendation& rec) {
  switch (rec.backend) {
    case LaplaceBackendKind::Diagonal:
      return std::make_unique<DiagonalBackend>();

    case LaplaceBackendKind::Tridiagonal:
      return std::make_unique<TridiagonalBackend>();

    case LaplaceBackendKind::Banded:
      return std::make_unique<BandedBackend>(rec.bandwidth);

    case LaplaceBackendKind::SparseLDLT:
      return std::make_unique<SparseLDLTBackend>();

    case LaplaceBackendKind::DenseLDLT:
      return std::make_unique<DenseLDLTBackend>();
  }

  throw std::runtime_error("Unknown Laplace backend recommendation");
}

inline std::unique_ptr<LaplaceBackend> CreateLaplaceBackendForHessian(
    const Eigen::SparseMatrix<double>& H,
    BackendRecommendation* out_recommendation = nullptr,
    const StructureDetectorOptions& options = StructureDetectorOptions()) {
  StructureDetector detector(options);
  BackendRecommendation rec = detector.Analyze(H);

  if (out_recommendation != nullptr) {
    *out_recommendation = rec;
  }

  auto backend = CreateLaplaceBackend(rec);
  backend->analyze_pattern(H);
  return backend;
}

}  // namespace laplace
}  // namespace quadra
EOF

cat > tests/test_laplace_backend_factory.cpp <<'EOF'
#include "../core/laplace/laplace_backend_factory.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <cmath>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <vector>

using quadra::laplace::BackendRecommendation;
using quadra::laplace::CreateLaplaceBackend;
using quadra::laplace::CreateLaplaceBackendForHessian;
using quadra::laplace::LaplaceBackendKind;
using quadra::laplace::LogDetDenseLDLT;
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
  H.makeCompressed();
  return H;
}

StructureDetectorOptions detector_options() {
  StructureDetectorOptions opts;
  opts.prefer_dense_for_small_matrices = false;
  opts.dense_size_cutoff = 0;
  opts.banded_width_cutoff = 8;
  opts.dense_fill_ratio = 0.40;
  return opts;
}

void test_backend_on_matrix(const Eigen::SparseMatrix<double>& H,
                            LaplaceBackendKind expected_backend,
                            const char* expected_name) {
  BackendRecommendation rec;
  auto backend = CreateLaplaceBackendForHessian(H, &rec, detector_options());

  expect_true(rec.backend == expected_backend, "backend recommendation mismatch");
  expect_true(std::string(backend->name()) == expected_name, "backend name mismatch");

  backend->factorize(H);

  expect_true(backend->is_spd(), "backend reports SPD");
  expect_near(backend->logdet(), LogDetDenseLDLT(Eigen::MatrixXd(H)), 1e-10,
              "backend logdet vs dense");
}

void test_diagonal_backend() {
  std::vector<Eigen::Triplet<double>> t;
  t.emplace_back(0, 0, 2.0);
  t.emplace_back(1, 1, 3.0);
  t.emplace_back(2, 2, 4.0);

  test_backend_on_matrix(make_sparse(3, t),
                         LaplaceBackendKind::Diagonal,
                         "diagonal");
}

void test_tridiagonal_backend() {
  const int n = 10;
  std::vector<Eigen::Triplet<double>> t;
  for (int i = 0; i < n; ++i) t.emplace_back(i, i, 4.0);
  for (int i = 0; i < n - 1; ++i) {
    t.emplace_back(i, i + 1, -0.25);
    t.emplace_back(i + 1, i, -0.25);
  }

  test_backend_on_matrix(make_sparse(n, t),
                         LaplaceBackendKind::Tridiagonal,
                         "tridiagonal");
}

void test_banded_backend() {
  const int n = 40;
  std::vector<Eigen::Triplet<double>> t;

  for (int i = 0; i < n; ++i) {
    t.emplace_back(i, i, 10.0);
    for (int d = 1; d <= 4; ++d) {
      if (i + d < n) {
        const double v = -0.05 / static_cast<double>(d);
        t.emplace_back(i, i + d, v);
        t.emplace_back(i + d, i, v);
      }
    }
  }

  test_backend_on_matrix(make_sparse(n, t),
                         LaplaceBackendKind::Banded,
                         "banded");
}

void test_sparse_backend() {
  const int n = 30;
  std::vector<Eigen::Triplet<double>> t;
  for (int i = 0; i < n; ++i) t.emplace_back(i, i, 5.0);

  t.emplace_back(0, 20, 0.1);
  t.emplace_back(20, 0, 0.1);
  t.emplace_back(3, 27, -0.1);
  t.emplace_back(27, 3, -0.1);

  test_backend_on_matrix(make_sparse(n, t),
                         LaplaceBackendKind::SparseLDLT,
                         "sparse_ldlt");
}

void test_dense_backend() {
  const int n = 8;
  std::vector<Eigen::Triplet<double>> t;

  Eigen::MatrixXd H = Eigen::MatrixXd::Constant(n, n, 0.05);
  H.diagonal().array() = 4.0;

  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      t.emplace_back(i, j, H(i, j));
    }
  }

  test_backend_on_matrix(make_sparse(n, t),
                         LaplaceBackendKind::DenseLDLT,
                         "dense_ldlt");
}

void test_sparse_symbolic_reuse() {
  const int n = 12;
  std::vector<Eigen::Triplet<double>> t1;
  std::vector<Eigen::Triplet<double>> t2;

  for (int i = 0; i < n; ++i) {
    t1.emplace_back(i, i, 5.0);
    t2.emplace_back(i, i, 6.0);
  }

  t1.emplace_back(0, 6, 0.1);
  t1.emplace_back(6, 0, 0.1);
  t2.emplace_back(0, 6, 0.2);
  t2.emplace_back(6, 0, 0.2);

  t1.emplace_back(2, 11, -0.1);
  t1.emplace_back(11, 2, -0.1);
  t2.emplace_back(2, 11, -0.2);
  t2.emplace_back(11, 2, -0.2);

  auto H1 = make_sparse(n, t1);
  auto H2 = make_sparse(n, t2);

  BackendRecommendation rec;
  auto backend = CreateLaplaceBackendForHessian(H1, &rec, detector_options());

  backend->factorize(H1);
  expect_true(backend->is_spd(), "symbolic reuse H1 SPD");
  const double logdet1 = backend->logdet();

  backend->factorize(H2);
  expect_true(backend->is_spd(), "symbolic reuse H2 SPD");
  const double logdet2 = backend->logdet();

  expect_near(logdet1, LogDetDenseLDLT(Eigen::MatrixXd(H1)), 1e-10,
              "reuse H1 logdet");
  expect_near(logdet2, LogDetDenseLDLT(Eigen::MatrixXd(H2)), 1e-10,
              "reuse H2 logdet");
}

}  // namespace

int main() {
  test_diagonal_backend();
  test_tridiagonal_backend();
  test_banded_backend();
  test_sparse_backend();
  test_dense_backend();
  test_sparse_symbolic_reuse();

  std::cout << "laplace backend factory tests passed\n";
  return 0;
}
EOF

cat > run_laplace_backend_factory_test.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/tests

set -x
"${CXX}" ${CXXFLAGS} \
  -Iexternal/Eigen \
  -I. \
  tests/test_laplace_backend_factory.cpp \
  -o build/tests/test_laplace_backend_factory

./build/tests/test_laplace_backend_factory
EOF

chmod +x run_laplace_backend_factory_test.sh

cat <<'EOF'

Installed Laplace backend interface + factory.

Run:
  ./run_laplace_backend_factory_test.sh

Next:
  wire BackendRecommendation + std::unique_ptr<LaplaceBackend>
  into PersistentLaplaceCache state, then adapt one example to use it.

EOF
