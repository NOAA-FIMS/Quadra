#include "../core/inference/ad_delta_method.hpp"
#include "../core/inference/delta_method.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <vector>

DECLARE_ADGRAPH()

int main() {
  const std::vector<double> theta{std::log(2.0), 0.35};

  Eigen::MatrixXd covariance(2, 2);
  covariance << 0.04, 0.01, 0.01, 0.09;

  auto fd_fun = [](const std::vector<double> &x) {
    const double a = std::exp(x[0]);
    const double p = 1.0 / (1.0 + std::exp(-x[1]));
    return a * p;
  };

  auto ad_fun = [](const std::vector<quadra::AD> &x) {
    const quadra::AD a = exp(x[0]);
    const quadra::AD p = 1.0 / (1.0 + exp(-x[1]));
    return a * p;
  };

  auto fd = quadra::delta_method_scalar(fd_fun, theta, covariance, 1.0e-6);

  auto ad = quadra::ad_delta_method_scalar(ad_fun, theta, covariance);

  if (!fd.success_m) {
    std::cerr << "FAIL: FD delta method failed: " << fd.message_m << "\n";
    return 1;
  }

  if (!ad.success_m) {
    std::cerr << "FAIL: AD delta method failed: " << ad.message_m << "\n";
    return 1;
  }

  const double estimate_error = std::abs(fd.estimate_m - ad.estimate_m);

  const double gradient_error =
      (fd.gradient_m - ad.gradient_m).cwiseAbs().maxCoeff();

  const double variance_error = std::abs(fd.variance_m - ad.variance_m);

  if (estimate_error > 1.0e-12) {
    std::cerr << "FAIL: estimate mismatch: " << estimate_error << "\n";
    return 1;
  }

  if (gradient_error > 1.0e-6) {
    std::cerr << "FAIL: gradient mismatch: " << gradient_error << "\n"
              << "FD: " << fd.gradient_m.transpose() << "\n"
              << "AD: " << ad.gradient_m.transpose() << "\n";
    return 1;
  }

  if (variance_error > 1.0e-6) {
    std::cerr << "FAIL: variance mismatch: " << variance_error << "\n"
              << "FD: " << fd.variance_m << "\n"
              << "AD: " << ad.variance_m << "\n";
    return 1;
  }

  std::cout << "PASS: AD delta method matches finite-difference delta method\n";
  std::cout << "  estimate: " << ad.estimate_m << "\n";
  std::cout << "  AD gradient: " << ad.gradient_m.transpose() << "\n";
  std::cout << "  FD gradient: " << fd.gradient_m.transpose() << "\n";
  std::cout << "  AD variance: " << ad.variance_m << "\n";
  std::cout << "  FD variance: " << fd.variance_m << "\n";

  return 0;
}
