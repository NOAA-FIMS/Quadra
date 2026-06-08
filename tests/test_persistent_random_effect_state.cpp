#include "../core/laplace/persistent_random_effect_state.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <stdexcept>

using quadra::laplace::NewtonSolveStatus;
using quadra::laplace::PersistentRandomEffectState;

void test_uninitialized_access_rejected() {
  PersistentRandomEffectState state;

  bool threw = false;
  try {
    (void)state.xhat();
  } catch (...) {
    threw = true;
  }

  if (!threw) {
    throw std::runtime_error("uninitialized xhat access was not rejected");
  }
}

void test_initialize_and_update() {
  PersistentRandomEffectState state;

  Eigen::VectorXd x = Eigen::VectorXd::Zero(3);
  x << 1.0, 2.0, 3.0;

  state.initialize(x);

  if (!state.initialized()) {
    throw std::runtime_error("state not initialized");
  }

  if (state.size() != 3) {
    throw std::runtime_error("wrong state size after initialize");
  }

  NewtonSolveStatus status;
  status.iterations = 4;
  status.objective = 12.5;
  status.grad_norm = 1e-8;
  status.converged = true;

  Eigen::VectorXd x2 = Eigen::VectorXd::Zero(3);
  x2 << 1.1, 2.1, 3.1;

  state.update_xhat(x2, status);

  if (state.status().iterations != 4 || !state.status().converged) {
    throw std::runtime_error("status not updated");
  }

  if (std::abs(state.xhat()[1] - 2.1) > 1e-12) {
    throw std::runtime_error("xhat not updated");
  }
}

void test_size_change_rejected_after_init() {
  PersistentRandomEffectState state;
  state.initialize(Eigen::VectorXd::Zero(3));

  NewtonSolveStatus status;

  bool threw = false;
  try {
    state.update_xhat(Eigen::VectorXd::Zero(4), status);
  } catch (...) {
    threw = true;
  }

  if (!threw) {
    throw std::runtime_error("size-changing update was not rejected");
  }
}

void test_resize_and_clear() {
  PersistentRandomEffectState state(5);

  if (state.initialized()) {
    throw std::runtime_error("resize constructor should not initialize xhat");
  }

  if (state.size() != 5) {
    throw std::runtime_error("wrong size after resize constructor");
  }

  state.initialize(Eigen::VectorXd::Ones(5));
  state.clear();

  if (state.initialized()) {
    throw std::runtime_error("state still initialized after clear");
  }

  if (state.size() != 0) {
    throw std::runtime_error("state size not cleared");
  }
}

int main() {
  test_uninitialized_access_rejected();
  test_initialize_and_update();
  test_size_change_rejected_after_init();
  test_resize_and_clear();

  std::cout << "persistent random effect state tests passed\\n";
  return 0;
}
