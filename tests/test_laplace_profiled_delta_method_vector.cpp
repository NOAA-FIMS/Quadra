#include "../core/laplace/laplace_profiled_delta_method_vector.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>

int main() {
  Eigen::VectorXd estimates(2);
  estimates << 10.0, 5.0;

  Eigen::MatrixXd g_theta(2, 2);
  g_theta << 1.0, -2.0, 0.5, 1.5;

  Eigen::MatrixXd g_u(2, 2);
  g_u << 3.0, 4.0, -1.0, 2.0;

  Eigen::MatrixXd du_dtheta(2, 2);
  du_dtheta << 0.5, -0.1, 0.2, 0.3;

  Eigen::MatrixXd theta_covariance(2, 2);
  theta_covariance << 0.04, 0.01, 0.01, 0.09;

  const auto result = quadra::compute_laplace_profiled_delta_method_vector(
      estimates, g_theta, g_u, du_dtheta, theta_covariance);

  if (!result.success_m) {
    std::cerr << "FAIL: vector profiled delta method failed: "
              << result.message_m << "\n";
    return 1;
  }

  const Eigen::MatrixXd expected_jacobian = g_theta + g_u * du_dtheta;

  const Eigen::MatrixXd expected_covariance =
      expected_jacobian * theta_covariance * expected_jacobian.transpose();

  const double jacobian_error =
      (result.jacobian_m - expected_jacobian).cwiseAbs().maxCoeff();

  const double covariance_error =
      (result.covariance_m - expected_covariance).cwiseAbs().maxCoeff();

  if (jacobian_error > 1.0e-12) {
    std::cerr << "FAIL: Jacobian mismatch\n";
    return 1;
  }

  if (covariance_error > 1.0e-12) {
    std::cerr << "FAIL: covariance mismatch\n";
    return 1;
  }

  if (result.std_error_m.size() != 2 || result.cv_m.size() != 2) {
    std::cerr << "FAIL: wrong output vector dimensions\n";
    return 1;
  }

  std::cout << "PASS: Laplace profiled vector delta-method utility\n";
  std::cout << "  Jacobian:\n" << result.jacobian_m << "\n";
  std::cout << "  covariance:\n" << result.covariance_m << "\n";
  std::cout << "  correlation:\n" << result.correlation_m << "\n";
  std::cout << "  std.errors: " << result.std_error_m.transpose() << "\n";
  std::cout << "  cvs: " << result.cv_m.transpose() << "\n";

  return 0;
}
