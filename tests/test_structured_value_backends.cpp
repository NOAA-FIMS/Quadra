#include "../core/laplace/structured_value_backend.hpp"
#include "../core/laplace/hessian_structure.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

using quadra::laplace::BandedValues;
using quadra::laplace::DiagonalValues;
using quadra::laplace::TridiagonalValues;
using quadra::laplace::LogDetSparseLDLT;
using quadra::laplace::logdet_banded_values_ldlt;
using quadra::laplace::logdet_diagonal_values;
using quadra::laplace::logdet_tridiagonal_values_ldlt;
using quadra::laplace::sparse_from_banded_values;
using quadra::laplace::sparse_from_diagonal_values;
using quadra::laplace::sparse_from_tridiagonal_values;

void expect_close(const double a, const double b, const char* msg) {
  const double diff = std::abs(a - b);
  const double scale = 1.0 + std::max(std::abs(a), std::abs(b));
  if (diff > 1e-9 * scale) {
    std::cerr << msg << ": a=" << a << " b=" << b
              << " diff=" << diff << "\n";
    throw std::runtime_error(msg);
  }
}

void test_diagonal() {
  DiagonalValues H;
  H.diag = Eigen::VectorXd(4);
  H.diag << 2.0, 3.0, 4.0, 5.0;

  const double direct = logdet_diagonal_values(H);
  const double sparse = LogDetSparseLDLT(sparse_from_diagonal_values(H));
  expect_close(direct, sparse, "diagonal logdet");
}

void test_tridiagonal() {
  TridiagonalValues H;
  H.diag = Eigen::VectorXd(5);
  H.offdiag = Eigen::VectorXd(4);

  H.diag << 4.0, 4.2, 4.4, 4.6, 4.8;
  H.offdiag << -0.5, -0.4, -0.3, -0.2;

  const double direct = logdet_tridiagonal_values_ldlt(H);
  const double sparse = LogDetSparseLDLT(sparse_from_tridiagonal_values(H));
  expect_close(direct, sparse, "tridiagonal logdet");
}

void test_banded() {
  BandedValues H;
  H.bandwidth = 2;
  H.diag = Eigen::VectorXd(6);
  H.diag << 5.0, 5.2, 5.4, 5.6, 5.8, 6.0;

  H.lower_bands.resize(2);
  H.lower_bands[0] = Eigen::VectorXd(5);
  H.lower_bands[1] = Eigen::VectorXd(4);
  H.lower_bands[0] << -0.4, -0.3, -0.25, -0.2, -0.1;
  H.lower_bands[1] << 0.08, 0.06, 0.04, 0.02;

  const double direct = logdet_banded_values_ldlt(H);
  const double sparse = LogDetSparseLDLT(sparse_from_banded_values(H));
  expect_close(direct, sparse, "banded logdet");
}

void test_non_spd_rejected() {
  DiagonalValues H;
  H.diag = Eigen::VectorXd(2);
  H.diag << 1.0, -1.0;

  bool threw = false;
  try {
    (void)logdet_diagonal_values(H);
  } catch (...) {
    threw = true;
  }

  if (!threw) {
    throw std::runtime_error("non-SPD diagonal was not rejected");
  }
}

int main() {
  test_diagonal();
  test_tridiagonal();
  test_banded();
  test_non_spd_rejected();

  std::cout << "structured value backend tests passed\n";
  return 0;
}
