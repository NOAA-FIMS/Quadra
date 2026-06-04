#include "../core/laplace/persistent_latent_state_runtime.hpp"
#include "../core/laplace/structured_value_backend.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <stdexcept>

using quadra::laplace::logdet_tridiagonal_values_ldlt;
using quadra::laplace::NewtonSolveStatus;
using quadra::laplace::PersistentLatentStateRuntime;
using quadra::laplace::TridiagonalValues;

void expect_close(const double a, const double b, const char *msg) {
  const double diff = std::abs(a - b);
  const double scale = 1.0 + std::max(std::abs(a), std::abs(b));

  if (diff > 1e-12 * scale) {
    std::cerr << msg << ": a=" << a << " b=" << b << " diff=" << diff << "\n";
    throw std::runtime_error(msg);
  }
}

TridiagonalValues make_tri(const int n) {
  TridiagonalValues H;
  H.diag = Eigen::VectorXd::Zero(n);
  H.offdiag = Eigen::VectorXd::Zero(std::max(0, n - 1));

  for (int i = 0; i < n; ++i) {
    H.diag[i] = 4.0 + 0.01 * i;
    if (i > 0) {
      H.offdiag[i - 1] = -0.2 + 0.001 * i;
    }
  }

  return H;
}

void test_empty_runtime() {
  PersistentLatentStateRuntime runtime;

  if (runtime.has_random_effects()) {
    throw std::runtime_error("new runtime unexpectedly has random effects");
  }

  if (runtime.has_structure()) {
    throw std::runtime_error("new runtime unexpectedly has structure");
  }
}

void test_random_effect_state_owned() {
  PersistentLatentStateRuntime runtime;

  Eigen::VectorXd x = Eigen::VectorXd::Zero(3);
  x << 1.0, 2.0, 3.0;

  NewtonSolveStatus status;
  status.iterations = 2;
  status.objective = 10.0;
  status.grad_norm = 1e-8;
  status.converged = true;

  runtime.random_effects().update_xhat(x, status);

  if (!runtime.has_random_effects()) {
    throw std::runtime_error("runtime did not report cached random effects");
  }

  if (runtime.random_effects().status().iterations != 2) {
    throw std::runtime_error("random effect status was not preserved");
  }
}

void test_structured_state_owned() {
  PersistentLatentStateRuntime runtime;

  const TridiagonalValues H = make_tri(20);
  runtime.structured().update_direct(H);

  if (!runtime.has_structure()) {
    throw std::runtime_error("runtime did not report cached structure");
  }

  expect_close(runtime.structured().logdet(), logdet_tridiagonal_values_ldlt(H),
               "runtime structured logdet");
}

void test_clear() {
  PersistentLatentStateRuntime runtime;

  NewtonSolveStatus status;
  runtime.random_effects().update_xhat(Eigen::VectorXd::Ones(4), status);
  runtime.structured().update_direct(make_tri(4));

  runtime.clear();

  if (runtime.has_random_effects()) {
    throw std::runtime_error("random effects not cleared");
  }

  if (runtime.has_structure()) {
    throw std::runtime_error("structure not cleared");
  }
}

int main() {
  test_empty_runtime();
  test_random_effect_state_owned();
  test_structured_state_owned();
  test_clear();

  std::cout << "persistent latent state runtime tests passed\n";
  return 0;
}
