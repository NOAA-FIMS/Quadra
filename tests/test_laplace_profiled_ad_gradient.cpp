#include "../core/laplace/laplace_profiled_ad_gradient.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <vector>

DECLARE_ADGRAPH()

int main() {
  const std::vector<double> fixed{2.0, -0.5};
  const std::vector<double> random{1.5, -2.0};

  auto derived = [](const std::vector<quadra::AD> &theta,
                    const std::vector<quadra::AD> &u) {
    return theta[0] * theta[0] + 3.0 * theta[1] + 2.0 * u[0] - 4.0 * u[1] +
           theta[0] * u[0];
  };

  const auto result =
      quadra::evaluate_profiled_ad_gradient_blocks(derived, fixed, random);

  if (!result.success_m) {
    std::cerr << "FAIL: profiled AD gradient failed: " << result.message_m
              << "\n";
    return 1;
  }

  const double expected_estimate = fixed[0] * fixed[0] + 3.0 * fixed[1] +
                                   2.0 * random[0] - 4.0 * random[1] +
                                   fixed[0] * random[0];

  Eigen::VectorXd expected_fixed(2);
  expected_fixed << 2.0 * fixed[0] + random[0], 3.0;

  Eigen::VectorXd expected_random(2);
  expected_random << 2.0 + fixed[0], -4.0;

  const double estimate_error = std::abs(result.estimate_m - expected_estimate);

  const double fixed_error = (result.gradient_fixed_m - expected_fixed).norm();

  const double random_error =
      (result.gradient_random_m - expected_random).norm();

  if (estimate_error > 1.0e-12 || fixed_error > 1.0e-12 ||
      random_error > 1.0e-12) {
    std::cerr << "FAIL: profiled AD gradient mismatch\n";
    std::cerr << "estimate error: " << estimate_error << "\n";
    std::cerr << "fixed gradient: " << result.gradient_fixed_m.transpose()
              << " expected " << expected_fixed.transpose() << "\n";
    std::cerr << "random gradient: " << result.gradient_random_m.transpose()
              << " expected " << expected_random.transpose() << "\n";
    return 1;
  }

  std::cout << "PASS: profiled AD gradient block test\n";
  std::cout << "  estimate: " << result.estimate_m << "\n";
  std::cout << "  g_theta: " << result.gradient_fixed_m.transpose() << "\n";
  std::cout << "  g_u: " << result.gradient_random_m.transpose() << "\n";

  return 0;
}
