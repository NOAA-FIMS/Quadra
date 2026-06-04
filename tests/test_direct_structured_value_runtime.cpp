#include "../core/laplace/persistent_structured_runtime.hpp"
#include "../core/laplace/structured_value_backend.hpp"

#include <cmath>
#include <iostream>
#include <stdexcept>

using quadra::laplace::BandedValues;
using quadra::laplace::DiagonalValues;
using quadra::laplace::LaplaceBackendKind;
using quadra::laplace::PersistentStructuredRuntimeState;
using quadra::laplace::TridiagonalValues;
using quadra::laplace::logdet_banded_values_ldlt;
using quadra::laplace::logdet_diagonal_values;
using quadra::laplace::logdet_tridiagonal_values_ldlt;

void expect_close(const double a, const double b, const char* msg) {
  const double diff = std::abs(a - b);
  const double scale = 1.0 + std::max(std::abs(a), std::abs(b));
  if (diff > 1e-9 * scale) {
    std::cerr << msg << ": a=" << a << " b=" << b
              << " diff=" << diff << "\n";
    throw std::runtime_error(msg);
  }
}

DiagonalValues make_diag(const int n, const double scale) {
  DiagonalValues H;
  H.diag = Eigen::VectorXd::Zero(n);
  for (int i = 0; i < n; ++i) {
    H.diag[i] = scale * (2.0 + 0.01 * i);
  }
  return H;
}

TridiagonalValues make_tri(const int n, const double scale) {
  TridiagonalValues H;
  H.diag = Eigen::VectorXd::Zero(n);
  H.offdiag = Eigen::VectorXd::Zero(std::max(0, n - 1));

  for (int i = 0; i < n; ++i) {
    H.diag[i] = scale * (4.0 + 0.01 * i);
    if (i > 0) {
      H.offdiag[i - 1] = scale * (-0.2 + 0.001 * i);
    }
  }
  return H;
}

BandedValues make_banded(const int n, const int bandwidth, const double scale) {
  BandedValues H;
  H.bandwidth = bandwidth;
  H.diag = Eigen::VectorXd::Zero(n);
  H.lower_bands.resize(static_cast<std::size_t>(bandwidth));

  for (int d = 1; d <= bandwidth; ++d) {
    H.lower_bands[static_cast<std::size_t>(d - 1)] =
        Eigen::VectorXd::Zero(std::max(0, n - d));
  }

  for (int i = 0; i < n; ++i) {
    double diag = scale * (10.0 + 0.01 * i);

    for (int d = 1; d <= bandwidth; ++d) {
      const int j = i - d;
      if (j < 0) continue;

      const double e =
          scale * (((d % 2 == 0) ? 0.015 : -0.025) /
                   static_cast<double>(d));

      H.lower_bands[static_cast<std::size_t>(d - 1)][j] = e;
      diag += 2.0 * std::abs(e);
    }

    H.diag[i] = diag;
  }

  return H;
}

void test_direct_diagonal() {
  PersistentStructuredRuntimeState state;
  const auto H = make_diag(20, 1.0);

  state.update_direct(H);

  if (state.backend_recommendation().backend != LaplaceBackendKind::Diagonal) {
    throw std::runtime_error("direct diagonal backend mismatch");
  }

  expect_close(state.logdet(), logdet_diagonal_values(H), "direct diagonal");
}

void test_direct_tridiagonal() {
  PersistentStructuredRuntimeState state;
  const auto H1 = make_tri(50, 1.0);
  const auto H2 = make_tri(50, 1.04);

  state.update_direct(H1);
  if (state.backend_recommendation().backend != LaplaceBackendKind::Tridiagonal) {
    throw std::runtime_error("direct tridiagonal backend mismatch");
  }
  expect_close(state.logdet(), logdet_tridiagonal_values_ldlt(H1),
               "direct tridiagonal first");

  state.update_direct(H2);
  expect_close(state.logdet(), logdet_tridiagonal_values_ldlt(H2),
               "direct tridiagonal second");
}

void test_direct_banded() {
  PersistentStructuredRuntimeState state;
  const auto H = make_banded(80, 5, 0.97);

  state.update_direct(H);

  if (state.backend_recommendation().backend != LaplaceBackendKind::Banded) {
    throw std::runtime_error("direct banded backend mismatch");
  }

  if (state.backend_recommendation().bandwidth != 5) {
    throw std::runtime_error("direct banded bandwidth mismatch");
  }

  expect_close(state.logdet(), logdet_banded_values_ldlt(H), "direct banded");
}

int main() {
  test_direct_diagonal();
  test_direct_tridiagonal();
  test_direct_banded();

  std::cout << "direct structured value runtime tests passed\n";
  return 0;
}
