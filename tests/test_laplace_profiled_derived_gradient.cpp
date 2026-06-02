#include "../core/laplace/laplace_profiled_derived_gradient.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>

int main() {
  Eigen::VectorXd g_theta(2);
  g_theta << 1.0, -2.0;

  Eigen::VectorXd g_u(2);
  g_u << 3.0, 4.0;

  Eigen::MatrixXd du_dtheta(2, 2);
  du_dtheta << 0.5, -0.1, 0.2, 0.3;

  const auto result = quadra::compute_laplace_profiled_derived_gradient(
      g_theta, g_u, du_dtheta);

  if (!result.success_m) {
    std::cerr << "FAIL: profiled derived gradient failed: " << result.message_m
              << "\n";
    return 1;
  }

  Eigen::VectorXd expected(2);
  expected = g_theta + (g_u.transpose() * du_dtheta).transpose();

  const double grad_error = (result.gradient_m - expected).norm();

  if (grad_error > 1.0e-12) {
    std::cerr << "FAIL: gradient mismatch\n";
    std::cerr << "actual: " << result.gradient_m.transpose() << "\n";
    std::cerr << "expected: " << expected.transpose() << "\n";
    return 1;
  }

  Eigen::MatrixXd covariance(2, 2);
  covariance << 0.04, 0.01, 0.01, 0.09;

  const double variance =
      quadra::delta_variance_from_gradient(result.gradient_m, covariance);

  if (!std::isfinite(variance) || variance <= 0.0) {
    std::cerr << "FAIL: bad propagated variance: " << variance << "\n";
    return 1;
  }

  std::cout << "PASS: Laplace profiled derived gradient utility\n";
  std::cout << "  gradient: " << result.gradient_m.transpose() << "\n";
  std::cout << "  variance: " << variance << "\n";
  std::cout << "  std.error: " << std::sqrt(variance) << "\n";

  return 0;
}
