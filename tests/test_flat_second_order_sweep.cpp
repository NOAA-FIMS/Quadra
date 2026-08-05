#include "../core/autodiff.hpp"

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <vector>

DECLARE_ADGRAPH();

namespace {

quadra::AD objective(const std::vector<quadra::AD> &x) {
  quadra::AD out = exp(x[0]) + x[0] * x[1];
  for (std::size_t i = 1; i < x.size(); ++i) {
    out += exp(0.1 * x[i]) + 0.3 * x[i - 1] * x[i] +
           0.05 * x[i] * x[i] * x[i];
  }
  return out;
}

Eigen::MatrixXd extract(const std::vector<quadra::AD> &x) {
  Eigen::MatrixXd hessian(x.size(), x.size());
  for (int i = 0; i < hessian.rows(); ++i) {
    for (int j = 0; j < hessian.cols(); ++j) {
      hessian(i, j) = had::GetAdjoint(x[i], x[j]);
    }
  }
  return hessian;
}

Eigen::MatrixXd fresh(const std::vector<double> &values) {
  quadra::TapeContext tape;
  quadra::ADScope scope(tape.graph);
  std::vector<quadra::AD> x = quadra::to_ad(values);
  quadra::AD y = objective(x);
  scope.backward(y);
  return extract(x);
}

} // namespace

int main() {
  std::vector<double> initial{0.2, -0.3, 0.4, 0.1, -0.2, 0.35};
  std::vector<double> updated{-0.1, 0.25, 0.3, -0.4, 0.15, 0.5};

  quadra::TapeContext tape;
  quadra::ADScope scope(tape.graph);
  std::vector<quadra::AD> x = quadra::to_ad(initial);
  quadra::AD y = objective(x);
  scope.backward(y);
  (void)tape.graph.FreezeHessianTopology();

  had::g_ADGraph = &tape.graph;
  for (std::size_t i = 0; i < x.size(); ++i) {
    had::SetValue(x[i], updated[i]);
  }
  tape.graph.Forward();
  tape.graph.ZeroAdjoints();
  had::SetAdjoint(y, 1.0);
  had::PropagateAdjoint();
  const Eigen::MatrixXd frozen = extract(x);
  const Eigen::MatrixXd reference = fresh(updated);

  const double error = (frozen - reference).cwiseAbs().maxCoeff();
  if (error > 1e-12) {
    std::cerr << "frozen-topology second-order sweep mismatch: " << error
              << "\n";
    return 1;
  }
  std::cout << "PASS: frozen-topology second-order sweep matches a fresh AD "
               "sweep\n";
  return 0;
}
