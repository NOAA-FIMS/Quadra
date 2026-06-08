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
make_sparse(int n, const std::vector<Eigen::Triplet<double>> &triplets) {
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

void test_backend_on_matrix(const Eigen::SparseMatrix<double> &H,
                            LaplaceBackendKind expected_backend,
                            const char *expected_name) {
  BackendRecommendation rec;
  auto backend = CreateLaplaceBackendForHessian(H, &rec, detector_options());

  expect_true(rec.backend == expected_backend,
              "backend recommendation mismatch");
  expect_true(std::string(backend->name()) == expected_name,
              "backend name mismatch");

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

  test_backend_on_matrix(make_sparse(3, t), LaplaceBackendKind::Diagonal,
                         "diagonal");
}

void test_tridiagonal_backend() {
  const int n = 10;
  std::vector<Eigen::Triplet<double>> t;
  for (int i = 0; i < n; ++i)
    t.emplace_back(i, i, 4.0);
  for (int i = 0; i < n - 1; ++i) {
    t.emplace_back(i, i + 1, -0.25);
    t.emplace_back(i + 1, i, -0.25);
  }

  test_backend_on_matrix(make_sparse(n, t), LaplaceBackendKind::Tridiagonal,
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

  test_backend_on_matrix(make_sparse(n, t), LaplaceBackendKind::Banded,
                         "banded");
}

void test_sparse_backend() {
  const int n = 30;
  std::vector<Eigen::Triplet<double>> t;
  for (int i = 0; i < n; ++i)
    t.emplace_back(i, i, 5.0);

  t.emplace_back(0, 20, 0.1);
  t.emplace_back(20, 0, 0.1);
  t.emplace_back(3, 27, -0.1);
  t.emplace_back(27, 3, -0.1);

  test_backend_on_matrix(make_sparse(n, t), LaplaceBackendKind::SparseLDLT,
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

  test_backend_on_matrix(make_sparse(n, t), LaplaceBackendKind::DenseLDLT,
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

} // namespace

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
