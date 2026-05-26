#include "../core/inference/ad_delta_method_vector.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <vector>

DECLARE_ADGRAPH()

int main() {
  const std::vector<double> theta{std::log(2.0), 0.35};

  Eigen::MatrixXd theta_covariance(2, 2);
  theta_covariance << 0.04, 0.01, 0.01, 0.09;

  auto f = [](const std::vector<quadra::AD> &x) {
    const quadra::AD y0 = exp(x[0]);
    const quadra::AD y1 = x[0] + 2.0 * x[1];
    return std::vector<quadra::AD>{y0, y1};
  };

  auto result = quadra::ad_delta_method_vector(f, theta, theta_covariance);

  if (!result.success_m) {
    std::cerr << "FAIL: AD vector delta method failed: " << result.message_m
              << "\n";
    return 1;
  }

  Eigen::VectorXd expected_estimate(2);
  expected_estimate << 2.0, theta[0] + 2.0 * theta[1];

  Eigen::MatrixXd expected_jacobian(2, 2);
  expected_jacobian << 2.0, 0.0, 1.0, 2.0;

  const Eigen::MatrixXd expected_covariance =
      expected_jacobian * theta_covariance * expected_jacobian.transpose();

  const double estimate_error =
      (result.estimate_m - expected_estimate).cwiseAbs().maxCoeff();

  const double jacobian_error =
      (result.jacobian_m - expected_jacobian).cwiseAbs().maxCoeff();

  const double covariance_error =
      (result.covariance_m - expected_covariance).cwiseAbs().maxCoeff();

  if (estimate_error > 1.0e-12) {
    std::cerr << "FAIL: estimate error too large: " << estimate_error << "\n";
    return 1;
  }

  if (jacobian_error > 1.0e-12) {
    std::cerr << "FAIL: Jacobian error too large: " << jacobian_error << "\n"
              << "estimated:\n"
              << result.jacobian_m << "\n"
              << "expected:\n"
              << expected_jacobian << "\n";
    return 1;
  }

  if (covariance_error > 1.0e-12) {
    std::cerr << "FAIL: covariance error too large: " << covariance_error
              << "\n"
              << "estimated:\n"
              << result.covariance_m << "\n"
              << "expected:\n"
              << expected_covariance << "\n";
    return 1;
  }

  std::cout << "PASS: AD vector delta-method test\n";
  std::cout << "  estimates: " << result.estimate_m.transpose() << "\n";
  std::cout << "  Jacobian:\n" << result.jacobian_m << "\n";
  std::cout << "  covariance:\n" << result.covariance_m << "\n";
  std::cout << "  std.errors: " << result.std_error_m.transpose() << "\n";

  return 0;
}
