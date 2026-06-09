#include "../core/uncertainty/reporting.hpp"

#include <Eigen/Core>

#include <cmath>
#include <iostream>
#include <vector>

int main() {
  Eigen::MatrixXd cov(3, 3);
  cov << 4.0, 1.0, 0.25,
         1.0, 9.0, 0.50,
         0.25, 0.50, 16.0;

  const auto corr = quadra::uncertainty::covariance_to_correlation_matrix(cov);

  if (std::abs(corr(0, 0) - 1.0) > 1.0e-12) return 1;
  if (std::abs(corr(0, 1) - (1.0 / 6.0)) > 1.0e-12) return 1;

  const auto diag = quadra::uncertainty::diagnose_covariance_matrix(cov);
  if (!diag.valid_covariance) {
    std::cerr << "expected valid covariance\n";
    return 1;
  }

  const auto decay = quadra::uncertainty::correlation_decay_summary(corr);
  if (decay.size() != 3) return 1;
  if (decay[0].count != 3) return 1;
  if (std::abs(decay[0].mean_correlation - 1.0) > 1.0e-12) return 1;

  Eigen::VectorXd log_est(2);
  log_est << std::log(10.0), std::log(20.0);

  Eigen::MatrixXd log_cov(2, 2);
  log_cov << 0.01, 0.002,
             0.002, 0.04;

  const auto natural_cov =
      quadra::uncertainty::lognormal_delta_covariance(log_est, log_cov);

  if (std::abs(natural_cov(0, 0) - 1.0) > 1.0e-12) return 1;
  if (std::abs(natural_cov(1, 1) - 16.0) > 1.0e-12) return 1;
  if (std::abs(natural_cov(0, 1) - 0.4) > 1.0e-12) return 1;

  std::vector<double> samples = {1.0, 2.0, 3.0, 4.0};
  auto row = quadra::uncertainty::summarize_samples(
      "scenario", 1, "quantity", 2.5, samples, "test");

  if (std::abs(row.mean - 2.5) > 1.0e-12) return 1;
  if (std::abs(row.median - 2.5) > 1.0e-12) return 1;
  if (row.n_samples != 4) return 1;

  std::cout << "uncertainty_reporting_test_passed\n";
  std::cout << "corr01=" << corr(0, 1) << "\n";
  std::cout << "min_eigenvalue=" << diag.min_eigenvalue << "\n";
  std::cout << "natural_cov01=" << natural_cov(0, 1) << "\n";

  return 0;
}
