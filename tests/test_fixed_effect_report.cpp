#include "../core/inference/fixed_effect_report.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

int main() {
  quadra::FixedEffectCovarianceResult cov;
  cov.success_m = true;
  cov.message_m = "ok";
  cov.covariance_m = Eigen::MatrixXd::Zero(2, 2);
  cov.covariance_m(0, 0) = 0.04;
  cov.covariance_m(1, 1) = 0.09;

  const std::vector<std::string> names{"theta_1", "theta_2"};
  const std::vector<double> estimates{1.0, -1.5};

  auto report = quadra::build_fixed_effect_report(names, estimates, cov);

  if (!report.success_m) {
    std::cerr << "FAIL: report failed: " << report.message_m << "\n";
    return 1;
  }

  if (report.rows_m.size() != 2) {
    std::cerr << "FAIL: wrong row count\n";
    return 1;
  }

  if (std::abs(report.rows_m[0].std_error_m - 0.2) > 1.0e-12) {
    std::cerr << "FAIL: wrong SE for theta_1\n";
    return 1;
  }

  if (std::abs(report.rows_m[1].std_error_m - 0.3) > 1.0e-12) {
    std::cerr << "FAIL: wrong SE for theta_2\n";
    return 1;
  }

  if (std::abs(report.rows_m[0].z_value_m - 5.0) > 1.0e-12) {
    std::cerr << "FAIL: wrong z for theta_1\n";
    return 1;
  }

  quadra::print_fixed_effect_report(report);

  std::cout << "PASS: fixed-effect report test\n";
  return 0;
}
