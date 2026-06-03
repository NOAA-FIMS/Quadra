#!/usr/bin/env bash
set -euo pipefail

# install_persistent_cache_backend_integration_v1.sh
#
# Integrates:
#   PersistentLaplaceCache
#   BackendRecommendation
#   LaplaceBackend
#
# into a reusable cache state scaffold.
#
# Adds:
#   core/laplace/persistent_laplace_runtime.hpp
#   tests/test_persistent_cache_backend_integration.cpp
#   run_persistent_cache_backend_integration_test.sh

mkdir -p core/laplace tests build/tests

cat > core/laplace/persistent_laplace_runtime.hpp <<'EOF'
#pragma once

#include "persistent_laplace_cache.hpp"
#include "laplace_backend_factory.hpp"

#include <memory>

namespace quadra {
namespace laplace {

struct RuntimeBackendState {
  bool backend_initialized = false;
  BackendRecommendation recommendation;
  std::unique_ptr<LaplaceBackend> backend;
};

template<class ModelAdapter>
class PersistentLaplaceRuntime {
 public:
  explicit PersistentLaplaceRuntime(ModelAdapter adapter)
      : cache_(std::move(adapter)) {}

  LaplaceEvaluation evaluate() {
    auto result = cache_.evaluate();

    if (!runtime_.backend_initialized) {
      Eigen::SparseMatrix<double> H =
          cache_.adapter().prototype_hessian();

      runtime_.backend =
          CreateLaplaceBackendForHessian(
              H,
              &runtime_.recommendation);

      runtime_.backend_initialized = true;
    }

    return result;
  }

  const RuntimeBackendState& runtime() const {
    return runtime_;
  }

  PersistentLaplaceCache<ModelAdapter>& cache() {
    return cache_;
  }

 private:
  PersistentLaplaceCache<ModelAdapter> cache_;
  RuntimeBackendState runtime_;
};

} // namespace laplace
} // namespace quadra
EOF

cat > tests/test_persistent_cache_backend_integration.cpp <<'EOF'
#include "../core/laplace/persistent_laplace_runtime.hpp"

#include <iostream>

using namespace quadra::laplace;

struct DummyAdapter {
  int random_size() const { return 5; }

  Eigen::VectorXd initial_random() const {
    return Eigen::VectorXd::Zero(5);
  }

  Eigen::VectorXd solve_random_effects(
      const Eigen::VectorXd& start) {
    return start;
  }

  LaplaceEvaluation evaluate_at_random(
      const Eigen::VectorXd&) {
    LaplaceEvaluation out;
    out.marginal = 1.0;
    return out;
  }

  Eigen::SparseMatrix<double> prototype_hessian() const {
    Eigen::SparseMatrix<double> H(5,5);
    H.insert(0,0)=1;
    H.insert(1,1)=1;
    H.insert(2,2)=1;
    H.insert(3,3)=1;
    H.insert(4,4)=1;
    return H;
  }
};

int main() {
  PersistentLaplaceRuntime<DummyAdapter> rt(
      DummyAdapter{});

  auto r = rt.evaluate();

  if (!rt.runtime().backend_initialized)
    throw std::runtime_error("backend not initialized");

  if (!rt.runtime().backend)
    throw std::runtime_error("backend missing");

  std::cout
      << "backend="
      << rt.runtime().backend->name()
      << "\n";

  std::cout
      << "persistent cache backend integration passed\n";

  return 0;
}
EOF

cat > run_persistent_cache_backend_integration_test.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/tests

c++ -std=c++17 -O2 -DNDEBUG -g \
  -Iexternal/Eigen \
  -I. \
  tests/test_persistent_cache_backend_integration.cpp \
  -o build/tests/test_persistent_cache_backend_integration

./build/tests/test_persistent_cache_backend_integration
EOF

chmod +x run_persistent_cache_backend_integration_test.sh

echo "Installed persistent cache/backend integration."
echo
echo "Run:"
echo "  ./run_persistent_cache_backend_integration_test.sh"
