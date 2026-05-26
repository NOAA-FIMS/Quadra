#include "../core/inference/delta_method.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <vector>

int main() {
  const std::vector<double> theta{1.0, 2.0};

  Eigen::MatrixXd covariance(2, 2);
  covariance << 0.04, 0.01, 0.01, 0.09;

  auto f = [](const std::vector<double> &x) { return x[0] + 2.0 * x[1]; };

  auto result = quadra::delta_method_scalar(f, theta, covariance, 1.0e-6);

  if (!result.success_m) {
    std::cerr << "FAIL: delta method failed: " << result.message_m << "\n";
    return 1;
  }

  const double expected_estimate = 5.0;
  const double expected_variance =
      1.0 * 1.0 * 0.04 + 2.0 * 2.0 * 0.09 + 2.0 * 1.0 * 2.0 * 0.01;

  if (std::abs(result.estimate_m - expected_estimate) > 1.0e-12) {
    std::cerr << "FAIL: wrong estimate: " << result.estimate_m << "\n";
    return 1;
  }

  if (std::abs(result.variance_m - expected_variance) > 1.0e-8) {
    std::cerr << "FAIL: wrong variance: " << result.variance_m << " expected "
              << expected_variance << "\n";
    return 1;
  }

  if (std::abs(result.gradient_m[0] - 1.0) > 1.0e-8 ||
      std::abs(result.gradient_m[1] - 2.0) > 1.0e-8) {
    std::cerr << "FAIL: wrong finite-difference gradient: "
              << result.gradient_m.transpose() << "\n";
    return 1;
  }

  std::cout << "PASS: scalar delta-method test\n";
  std::cout << "  estimate: " << result.estimate_m << "\n";
  std::cout << "  variance: " << result.variance_m << "\n";
  std::cout << "  std.error: " << result.std_error_m << "\n";

  return 0;
}
