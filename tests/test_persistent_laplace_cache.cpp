#include "../core/laplace/persistent_laplace_cache.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <stdexcept>

using quadra::laplace::LaplaceEvaluation;
using quadra::laplace::PersistentLaplaceCache;

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

// Toy adapter:
//   joint(u) = 0.5 * sum_i h_i * (u_i - mu_i)^2
//
// The optimum is uhat = mu.
// H is diagonal with logdet = sum log(h_i).
struct ToyDiagonalLaplaceAdapter {
  Eigen::VectorXd mu;
  Eigen::VectorXd h;
  int solve_calls = 0;
  int eval_calls = 0;

  int random_size() const { return static_cast<int>(mu.size()); }

  Eigen::VectorXd initial_random() const {
    return Eigen::VectorXd::Zero(mu.size());
  }

  Eigen::VectorXd solve_random_effects(const Eigen::VectorXd &start) {
    ++solve_calls;

    // One Newton step for this quadratic.
    Eigen::VectorXd grad = h.array() * (start - mu).array();
    Eigen::VectorXd step = grad.array() / h.array();

    return start - step;
  }

  LaplaceEvaluation evaluate_at_random(const Eigen::VectorXd &uhat) {
    ++eval_calls;

    const Eigen::VectorXd residual = uhat - mu;

    LaplaceEvaluation out;
    out.random_size = static_cast<int>(mu.size());
    out.joint = 0.5 * (h.array() * residual.array().square()).sum();
    out.logdet = h.array().log().sum();
    out.correction =
        0.5 * out.logdet -
        0.5 * static_cast<double>(out.random_size) * std::log(2.0 * M_PI);
    out.marginal = out.joint + out.correction;
    out.grad_norm = (h.array() * residual.array()).matrix().norm();
    out.nnz_hessian = out.random_size;
    out.backend = "diagonal";
    return out;
  }
};

} // namespace

int main() {
  ToyDiagonalLaplaceAdapter adapter;
  adapter.mu = Eigen::VectorXd::LinSpaced(5, -1.0, 1.0);
  adapter.h = Eigen::VectorXd::LinSpaced(5, 2.0, 6.0);

  PersistentLaplaceCache<ToyDiagonalLaplaceAdapter> cache(adapter);

  expect_true(!cache.state().initialized, "cache starts uninitialized");

  const auto cold = cache.evaluate_cold();
  expect_true(cache.state().initialized, "cache initialized after cold");
  expect_true(cache.state().cold_evaluations == 1, "cold eval count");
  expect_true(cache.state().warm_evaluations == 0, "warm eval count");
  expect_true(cache.state().cached_evaluations == 0, "cached eval count");
  expect_near(cold.grad_norm, 0.0, 1e-12, "cold grad norm");
  expect_near((cache.uhat() - adapter.mu).norm(), 0.0, 1e-12,
              "cached uhat equals optimum");

  const auto cached = cache.evaluate_cached_no_solve();
  expect_true(cache.state().cached_evaluations == 1,
              "cached eval count after cached");
  expect_near(cached.marginal, cold.marginal, 1e-12,
              "cached marginal equals cold");

  const auto warm = cache.evaluate_warm();
  expect_true(cache.state().warm_evaluations == 1,
              "warm eval count after warm");
  expect_near(warm.marginal, cold.marginal, 1e-12, "warm marginal equals cold");
  expect_near(warm.grad_norm, 0.0, 1e-12, "warm grad norm");

  cache.reset();
  expect_true(!cache.state().initialized, "cache reset uninitialized");

  const auto eval_default_1 = cache.evaluate();
  const auto eval_default_2 = cache.evaluate();
  expect_true(cache.state().cold_evaluations == 2,
              "evaluate first call uses cold");
  expect_true(cache.state().warm_evaluations == 2,
              "evaluate second call uses warm");
  expect_near(eval_default_1.marginal, eval_default_2.marginal, 1e-12,
              "default evaluations agree");

  std::cout << "persistent Laplace cache tests passed\n";
  return 0;
}
