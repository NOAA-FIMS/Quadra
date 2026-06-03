#!/usr/bin/env bash
set -euo pipefail

# install_persistent_laplace_cache_scaffold_v1.sh
#
# Adds a reusable, model-agnostic persistent Laplace cache scaffold.
#
# This is intentionally lightweight and non-invasive:
#
#   core/laplace/persistent_laplace_cache.hpp
#
# It provides:
#   - LaplaceCacheState
#   - LaplaceEvaluation
#   - PersistentLaplaceCache<ModelAdapter>
#
# ModelAdapter must provide:
#   int random_size() const;
#   Eigen::VectorXd initial_random() const;
#   Eigen::VectorXd solve_random_effects(const Eigen::VectorXd& start);
#   LaplaceEvaluation evaluate_at_random(const Eigen::VectorXd& uhat);
#
# The adapter owns all model-specific details:
#   - objective
#   - Hessian construction
#   - structure-aware logdet
#   - backend choice
#
# The cache owns:
#   - previous uhat
#   - initialized flag
#   - cold/warm/cached semantics
#
# Adds:
#   tests/test_persistent_laplace_cache.cpp
#   run_persistent_laplace_cache_test.sh
#
# This gives us the core architecture target without touching existing examples.

mkdir -p core/laplace tests build/tests

cat > core/laplace/persistent_laplace_cache.hpp <<'EOF'
#pragma once

#include <Eigen/Dense>

#include <cmath>
#include <stdexcept>
#include <string>

namespace quadra {
namespace laplace {

struct LaplaceEvaluation {
  double marginal = 0.0;
  double joint = 0.0;
  double logdet = 0.0;
  double correction = 0.0;
  double grad_norm = 0.0;
  int random_size = 0;
  int nnz_hessian = 0;
  std::string backend = "unknown";
};

struct LaplaceCacheState {
  bool initialized = false;
  Eigen::VectorXd uhat;
  int cold_evaluations = 0;
  int warm_evaluations = 0;
  int cached_evaluations = 0;
};

template <class ModelAdapter>
class PersistentLaplaceCache {
 public:
  explicit PersistentLaplaceCache(ModelAdapter adapter)
      : adapter_(std::move(adapter)) {
    state_.uhat = adapter_.initial_random();
    if (state_.uhat.size() != adapter_.random_size()) {
      throw std::invalid_argument(
          "PersistentLaplaceCache: initial_random size does not match random_size");
    }
  }

  const LaplaceCacheState& state() const { return state_; }

  const Eigen::VectorXd& uhat() const {
    if (!state_.initialized) {
      throw std::runtime_error("PersistentLaplaceCache: uhat requested before initialization");
    }
    return state_.uhat;
  }

  void reset() {
    state_.initialized = false;
    state_.uhat = adapter_.initial_random();
  }

  // Full cold evaluation from adapter-provided initial random effects.
  LaplaceEvaluation evaluate_cold() {
    state_.uhat = adapter_.solve_random_effects(adapter_.initial_random());
    state_.initialized = true;
    ++state_.cold_evaluations;
    return adapter_.evaluate_at_random(state_.uhat);
  }

  // Warm-start solve from the previously cached uhat.
  LaplaceEvaluation evaluate_warm() {
    if (!state_.initialized) {
      return evaluate_cold();
    }

    state_.uhat = adapter_.solve_random_effects(state_.uhat);
    ++state_.warm_evaluations;
    return adapter_.evaluate_at_random(state_.uhat);
  }

  // Evaluate Laplace at cached uhat without re-solving.
  LaplaceEvaluation evaluate_cached_no_solve() {
    if (!state_.initialized) {
      throw std::runtime_error(
          "PersistentLaplaceCache: cached evaluation requested before initialization");
    }

    ++state_.cached_evaluations;
    return adapter_.evaluate_at_random(state_.uhat);
  }

  // Common production default:
  //   first call cold-solves,
  //   subsequent calls warm-start.
  LaplaceEvaluation evaluate() {
    if (!state_.initialized) {
      return evaluate_cold();
    }
    return evaluate_warm();
  }

  ModelAdapter& adapter() { return adapter_; }
  const ModelAdapter& adapter() const { return adapter_; }

 private:
  ModelAdapter adapter_;
  LaplaceCacheState state_;
};

}  // namespace laplace
}  // namespace quadra
EOF

cat > tests/test_persistent_laplace_cache.cpp <<'EOF'
#include "../core/laplace/persistent_laplace_cache.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <stdexcept>

using quadra::laplace::LaplaceEvaluation;
using quadra::laplace::PersistentLaplaceCache;

namespace {

void expect_true(bool cond, const char* msg) {
  if (!cond) {
    throw std::runtime_error(msg);
  }
}

void expect_near(double a, double b, double tol, const char* msg) {
  if (std::abs(a - b) > tol) {
    std::cerr << msg << ": a=" << a << " b=" << b
              << " diff=" << std::abs(a - b) << "\n";
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

  int random_size() const {
    return static_cast<int>(mu.size());
  }

  Eigen::VectorXd initial_random() const {
    return Eigen::VectorXd::Zero(mu.size());
  }

  Eigen::VectorXd solve_random_effects(const Eigen::VectorXd& start) {
    ++solve_calls;

    // One Newton step for this quadratic.
    Eigen::VectorXd grad = h.array() * (start - mu).array();
    Eigen::VectorXd step = grad.array() / h.array();

    return start - step;
  }

  LaplaceEvaluation evaluate_at_random(const Eigen::VectorXd& uhat) {
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

}  // namespace

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
  expect_near((cache.uhat() - adapter.mu).norm(), 0.0, 1e-12, "cached uhat equals optimum");

  const auto cached = cache.evaluate_cached_no_solve();
  expect_true(cache.state().cached_evaluations == 1, "cached eval count after cached");
  expect_near(cached.marginal, cold.marginal, 1e-12, "cached marginal equals cold");

  const auto warm = cache.evaluate_warm();
  expect_true(cache.state().warm_evaluations == 1, "warm eval count after warm");
  expect_near(warm.marginal, cold.marginal, 1e-12, "warm marginal equals cold");
  expect_near(warm.grad_norm, 0.0, 1e-12, "warm grad norm");

  cache.reset();
  expect_true(!cache.state().initialized, "cache reset uninitialized");

  const auto eval_default_1 = cache.evaluate();
  const auto eval_default_2 = cache.evaluate();
  expect_true(cache.state().cold_evaluations == 2, "evaluate first call uses cold");
  expect_true(cache.state().warm_evaluations == 2, "evaluate second call uses warm");
  expect_near(eval_default_1.marginal, eval_default_2.marginal, 1e-12,
              "default evaluations agree");

  std::cout << "persistent Laplace cache tests passed\n";
  return 0;
}
EOF

cat > run_persistent_laplace_cache_test.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/tests

set -x
"${CXX}" ${CXXFLAGS} \
  -Iexternal/Eigen \
  -I. \
  tests/test_persistent_laplace_cache.cpp \
  -o build/tests/test_persistent_laplace_cache

./build/tests/test_persistent_laplace_cache
EOF

chmod +x run_persistent_laplace_cache_test.sh

cat <<'EOF'

Installed persistent Laplace cache scaffold.

Run:
  ./run_persistent_laplace_cache_test.sh

Next:
  adapt the age-structured flat-band example to use PersistentLaplaceCache<ModelAdapter>,
  then wire the same pattern into the real Laplace evaluator.

EOF
