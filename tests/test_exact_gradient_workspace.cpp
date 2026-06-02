#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../core/laplace/exact_gradient_workspace.hpp"

DECLARE_ADGRAPH()

namespace {

had::AReal objective(const std::vector<had::AReal> &theta,
                     const std::vector<had::AReal> &u) {
  const had::AReal &a = theta[0];
  const had::AReal &b = theta[1];

  had::AReal f(0.0);

  for (std::size_t i = 0; i < u.size(); ++i) {
    f = f + 0.5 * (u[i] - a) * (u[i] - a) + exp(b) * u[i] * u[i];
  }

  for (std::size_t i = 1; i < u.size(); ++i) {
    const had::AReal diff = u[i] - u[i - 1];
    f = f + 0.25 * diff * diff;
  }

  return f;
}

} // namespace

int main() {
  using quadra::laplace::ExactGradientWorkspace;
  using quadra::laplace::MakeTridiagonalHdotPattern;

  std::vector<had::AReal> theta(2);
  std::vector<had::AReal> u(4);

  ExactGradientWorkspace workspace;

  had::AReal f = workspace.Build(
      [&]() {
        theta[0] = had::AReal(0.2);
        theta[1] = had::AReal(-0.5);

        for (std::size_t i = 0; i < u.size(); ++i) {
          u[i] = had::AReal(0.1 * static_cast<double>(i + 1));
        }

        return objective(theta, u);
      },
      &theta, &u);

  (void)f;

  workspace.PropagateBaseAdjoint();

  workspace.SeedTotalDirections(2, [](std::size_t k,
                                      Eigen::VectorXd &theta_direction,
                                      Eigen::VectorXd &random_direction) {
    theta_direction = Eigen::VectorXd::Zero(2);
    random_direction = Eigen::VectorXd::Zero(4);

    theta_direction[static_cast<int>(k)] = 1.0;

    for (int i = 0; i < random_direction.size(); ++i) {
      random_direction[i] = 0.01 * static_cast<double>((i + 1) * (k + 1));
    }
  });

  workspace.PropagateDirectionalBatch();

  const auto pattern = MakeTridiagonalHdotPattern(4);

  Eigen::MatrixXd Hinv = Eigen::MatrixXd::Identity(4, 4);

  const Eigen::VectorXd traces = workspace.TraceTerms(Hinv, pattern);

  const Eigen::VectorXd selected_inverse_traces =
      workspace.TraceTermsSelectedInverse(
          [&](int row, int col) { return Hinv(row, col); }, pattern);

  if (traces.size() != 2) {
    throw std::runtime_error("TraceTerms returned wrong number of directions");
  }

  if (selected_inverse_traces.size() != traces.size()) {
    throw std::runtime_error(
        "TraceTermsSelectedInverse returned wrong number of directions");
  }

  if ((traces - selected_inverse_traces).cwiseAbs().maxCoeff() > 1.0e-12) {
    throw std::runtime_error(
        "TraceTermsSelectedInverse does not match dense TraceTerms");
  }

  if (!std::isfinite(traces[0]) || !std::isfinite(traces[1])) {
    throw std::runtime_error("TraceTerms returned non-finite values");
  }

  const Eigen::VectorXd joint_gradient = Eigen::VectorXd::Zero(2);

  const auto assembled = workspace.AssembleExactGradient(
      12.5, 1.25, joint_gradient,
      [&](int row, int col) { return Hinv(row, col); }, pattern);

  if (std::abs(assembled.objective - 13.125) > 1.0e-12) {
    throw std::runtime_error("AssembleExactGradient returned wrong objective");
  }

  if (assembled.gradient.size() != 2 || assembled.trace_terms.size() != 2) {
    throw std::runtime_error(
        "AssembleExactGradient returned wrong vector sizes");
  }

  if ((assembled.gradient - 0.5 * traces).cwiseAbs().maxCoeff() > 1.0e-12) {
    throw std::runtime_error(
        "AssembleExactGradient gradient does not match trace assembly");
  }

  const Eigen::MatrixXd hdot0 = workspace.ExtractHdotDense(0, pattern);

  if (hdot0.rows() != 4 || hdot0.cols() != 4) {
    throw std::runtime_error("ExtractHdotDense returned wrong dimensions");
  }

  const auto triplets = workspace.ExtractHdotTriplets(1, pattern);

  if (triplets.empty()) {
    throw std::runtime_error("ExtractHdotTriplets returned empty result");
  }

  std::cout << "exact gradient workspace tests passed\n";
  std::cout << "vertices = " << workspace.HadWorkspace().VertexCount() << "\n";
  std::cout << "trace 0 = " << traces[0] << "\n";
  std::cout << "trace 1 = " << traces[1] << "\n";
  std::cout << "triplets = " << triplets.size() << "\n";

  return 0;
}
