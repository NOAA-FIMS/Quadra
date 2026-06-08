#include "../core/laplace/hessian_structure.hpp"
#include "../core/laplace/laplace_backend.hpp"

#include <Eigen/Sparse>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

using quadra::laplace::DiagonalBackend;
using quadra::laplace::LogDetSparseLDLT;
using quadra::laplace::TridiagonalBackend;

void expect_close(const double a, const double b, const char *msg) {
  const double diff = std::abs(a - b);
  const double scale = 1.0 + std::max(std::abs(a), std::abs(b));

  if (diff > 1e-9 * scale) {
    std::cerr << msg << ": a=" << a << " b=" << b << " diff=" << diff << "\n";
    throw std::runtime_error(msg);
  }
}

Eigen::SparseMatrix<double> make_diagonal() {
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.emplace_back(0, 0, 2.0);
  triplets.emplace_back(1, 1, 3.0);
  triplets.emplace_back(2, 2, 4.0);
  triplets.emplace_back(3, 3, 5.0);

  Eigen::SparseMatrix<double> H(4, 4);
  H.setFromTriplets(triplets.begin(), triplets.end());
  H.makeCompressed();
  return H;
}

Eigen::SparseMatrix<double> make_tridiagonal() {
  const int n = 5;
  std::vector<Eigen::Triplet<double>> triplets;

  for (int i = 0; i < n; ++i) {
    triplets.emplace_back(i, i, 4.0 + 0.2 * i);

    if (i > 0) {
      const double e = -0.1 * i;
      triplets.emplace_back(i, i - 1, e);
      triplets.emplace_back(i - 1, i, e);
    }
  }

  Eigen::SparseMatrix<double> H(n, n);
  H.setFromTriplets(triplets.begin(), triplets.end());
  H.makeCompressed();
  return H;
}

void test_diagonal_backend_uses_value_path() {
  const Eigen::SparseMatrix<double> H = make_diagonal();

  DiagonalBackend backend;
  backend.factorize(H);

  if (!backend.is_spd()) {
    throw std::runtime_error("diagonal backend reported non-SPD");
  }

  expect_close(backend.logdet(), LogDetSparseLDLT(H),
               "diagonal backend logdet");
}

void test_tridiagonal_backend_uses_value_path() {
  const Eigen::SparseMatrix<double> H = make_tridiagonal();

  TridiagonalBackend backend;
  backend.factorize(H);

  if (!backend.is_spd()) {
    throw std::runtime_error("tridiagonal backend reported non-SPD");
  }

  expect_close(backend.logdet(), LogDetSparseLDLT(H),
               "tridiagonal backend logdet");
}

int main() {
  test_diagonal_backend_uses_value_path();
  test_tridiagonal_backend_uses_value_path();

  std::cout << "laplace backend structured-value tests passed\n";
  return 0;
}
