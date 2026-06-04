#include "../core/laplace/persistent_latent_state_runtime.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <stdexcept>

using quadra::laplace::NewtonSolveStatus;
using quadra::laplace::PersistentLatentStateRuntime;

struct ToyContext {
  Eigen::VectorXd target;
};

struct ToySolveResult {
  Eigen::VectorXd xhat;
  NewtonSolveStatus status;
};

struct ToySolverPolicy {
  int initial_calls = 0;
  int solve_calls = 0;
  Eigen::VectorXd last_initial_x;

  Eigen::VectorXd initial_x(ToyContext& ctx) {
    ++initial_calls;
    return Eigen::VectorXd::Zero(ctx.target.size());
  }

  ToySolveResult solve(ToyContext& ctx, const Eigen::VectorXd& initial_x) {
    ++solve_calls;
    last_initial_x = initial_x;

    ToySolveResult out;
    out.xhat = ctx.target;
    out.status.iterations = 1;
    out.status.objective = (ctx.target - initial_x).squaredNorm();
    out.status.grad_norm = 0.0;
    out.status.converged = true;
    return out;
  }
};

void test_solver_policy_initializes_then_warm_starts() {
  PersistentLatentStateRuntime runtime;
  ToyContext ctx;
  ctx.target = Eigen::VectorXd::Zero(3);
  ctx.target << 1.0, 2.0, 3.0;

  ToySolverPolicy solver;

  const auto first = runtime.solve_random_effects(ctx, solver);

  if (!runtime.has_random_effects()) {
    throw std::runtime_error("runtime did not cache random effects");
  }

  if (solver.initial_calls != 1 || solver.solve_calls != 1) {
    throw std::runtime_error("first solve did not call initial/solve exactly once");
  }

  if ((runtime.random_effects().xhat() - ctx.target).norm() > 1e-12) {
    throw std::runtime_error("runtime did not store first xhat");
  }

  ctx.target << 1.5, 2.5, 3.5;
  const auto second = runtime.solve_random_effects(ctx, solver);

  if (solver.initial_calls != 1) {
    throw std::runtime_error("second solve should not call initial_x");
  }

  if (solver.solve_calls != 2) {
    throw std::runtime_error("second solve not called");
  }

  if ((solver.last_initial_x - first.xhat).norm() > 1e-12) {
    throw std::runtime_error("second solve was not warm-started from cached xhat");
  }

  if ((runtime.random_effects().xhat() - ctx.target).norm() > 1e-12) {
    throw std::runtime_error("runtime did not store second xhat");
  }

  (void)second;
}

int main() {
  test_solver_policy_initializes_then_warm_starts();

  std::cout << "persistent latent solver policy tests passed\n";
  return 0;
}
