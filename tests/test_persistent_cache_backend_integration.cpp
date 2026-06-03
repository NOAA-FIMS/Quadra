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
