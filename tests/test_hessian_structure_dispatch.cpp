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

void expect_true(bool cond, const char *msg) {
  if (!cond) {
    throw std::runtime_error(msg);
  }
}

void expect_near(double a, double b, double tol, const char *msg) {
  if (std::abs(a - b) > tol) {
    std::cerr << msg << ": a=" << a << " b=" << b << " diff=" << std::abs(a - b)
              << "\n";
    throw std::runtime_error(msg);
  }
}

Eigen::SparseMatrix<double>
make_sparse_from_triplets(int n,
                          const std::vector<Eigen::Triplet<double>> &triplets) {
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
  for (int i = 0; i < n; ++i)
    t.emplace_back(i, i, 4.0);
  for (int i = 0; i < n - 1; ++i) {
    t.emplace_back(i, i + 1, -1.0);
    t.emplace_back(i + 1, i, -1.0);
  }

  auto H = make_sparse_from_triplets(n, t);
  StructureInfo info;
  HessianStructure backend;
  const double auto_logdet =
      AutomaticLogDet(H, StructureOptions(), &backend, &info);

  const double dense_logdet =
      quadra::laplace::LogDetDenseLDLT(Eigen::MatrixXd(H));

  expect_true(info.detected == HessianStructure::Tridiagonal,
              "tridiagonal detected");
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
  for (int i = 0; i < n; ++i)
    t.emplace_back(i, i, 5.0);

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

  expect_true(info.detected == HessianStructure::SparsePattern,
              "sparse detected");
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

} // namespace

int main() {
  test_diagonal();
  test_tridiagonal();
  test_banded();
  test_sparse_pattern();
  test_dense();

  std::cout << "hessian structure dispatch tests passed\n";
  return 0;
}
