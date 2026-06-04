#include "../core/laplace/laplace_evaluator.hpp"
#include "../core/laplace/structured_value_backend.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <stdexcept>

using quadra::laplace::BandedValues;
using quadra::laplace::LaplaceEvaluator;
using quadra::laplace::NewtonSolveStatus;
using quadra::laplace::PersistentStructuredRuntimeState;
using quadra::laplace::TridiagonalValues;

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

  Eigen::VectorXd initial_x(ToyContext& ctx) {
    ++initial_calls;
    return Eigen::VectorXd::Zero(ctx.target.size());
  }

  ToySolveResult solve(ToyContext& ctx, const Eigen::VectorXd& initial_x) {
    ++solve_calls;

    ToySolveResult out;
    out.xhat = ctx.target;
    out.status.iterations = 1;
    out.status.objective = (ctx.target - initial_x).squaredNorm();
    out.status.grad_norm = 0.0;
    out.status.converged = true;
    return out;
  }
};

struct ToyEvalResult {
  double objective = 0.0;
  double logdet = 0.0;
};

struct ToyEvaluationPolicy {
  ToyEvalResult evaluate(ToyContext& ctx,
                         const Eigen::VectorXd& xhat,
                         PersistentStructuredRuntimeState& structured) {
    TridiagonalValues H;
    const int n = static_cast<int>(xhat.size());
    H.diag = Eigen::VectorXd::Constant(n, 4.0);
    H.offdiag = Eigen::VectorXd::Constant(std::max(0, n - 1), -0.2);

    structured.update_direct(H);

    ToyEvalResult out;
    out.logdet = structured.logdet();
    out.objective = (xhat - ctx.target).squaredNorm() + out.logdet;
    return out;
  }
};

void test_laplace_evaluator_owns_runtime() {
  ToyContext ctx;
  ctx.target = Eigen::VectorXd::Zero(3);
  ctx.target << 1.0, 2.0, 3.0;

  LaplaceEvaluator<ToyContext, ToySolverPolicy, ToyEvaluationPolicy> evaluator;

  const ToyEvalResult first = evaluator.evaluate(ctx);

  if (!evaluator.runtime().has_random_effects()) {
    throw std::runtime_error("evaluator did not cache random effects");
  }

  if (!evaluator.runtime().has_structure()) {
    throw std::runtime_error("evaluator did not cache structure");
  }

  if (!(first.logdet > 0.0)) {
    throw std::runtime_error("unexpected first logdet");
  }

  ctx.target << 1.5, 2.5, 3.5;

  const ToyEvalResult second = evaluator.evaluate(ctx);

  if (!(second.logdet > 0.0)) {
    throw std::runtime_error("unexpected second logdet");
  }

  if (evaluator.solver_policy().initial_calls != 1) {
    throw std::runtime_error("solver initial_x should only be called once");
  }

  if (evaluator.solver_policy().solve_calls != 2) {
    throw std::runtime_error("solver solve call count mismatch");
  }
}

void test_clear_resets_runtime() {
  ToyContext ctx;
  ctx.target = Eigen::VectorXd::Ones(2);

  LaplaceEvaluator<ToyContext, ToySolverPolicy, ToyEvaluationPolicy> evaluator;
  (void)evaluator.evaluate(ctx);

  evaluator.clear();

  if (evaluator.runtime().has_random_effects()) {
    throw std::runtime_error("clear did not reset random effects");
  }

  if (evaluator.runtime().has_structure()) {
    throw std::runtime_error("clear did not reset structure");
  }
}

int main() {
  test_laplace_evaluator_owns_runtime();
  test_clear_resets_runtime();

  std::cout << "laplace evaluator tests passed\n";
  return 0;
}
