#include "../core/laplace/functional_analysis_report.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <stdexcept>

void expect_close(double actual, double expected, const char *message) {
  if (std::abs(actual - expected) > 1.0e-12) {
    throw std::runtime_error(message);
  }
}

int main() {
  Eigen::MatrixXd hessian = Eigen::MatrixXd::Zero(2, 2);
  hessian(0, 0) = 4.0;
  hessian(1, 1) = 1.0;
  const auto spectrum = quadra::summarize_spectral_structure(hessian);

  if (!spectrum.available || spectrum.eigenvalues_desc.size() != 2 ||
      spectrum.signed_eigenvalues_desc.size() != 2 ||
      spectrum.cumulative_share.size() != 2 ||
      spectrum.positive_eigen_count != 2 ||
      spectrum.negative_eigen_count != 0) {
    throw std::runtime_error("full spectral payload unavailable");
  }
  expect_close(spectrum.eigenvalues_desc[0], 4.0,
               "largest eigenvalue mismatch");
  expect_close(spectrum.eigenvalues_desc[1], 1.0,
               "smallest eigenvalue mismatch");
  expect_close(spectrum.cumulative_share[0], 0.8,
               "cumulative share mismatch");
  expect_close(spectrum.participation_ratio, 25.0 / 17.0,
               "participation ratio mismatch");
  expect_close(spectrum.stable_rank, 17.0 / 16.0,
               "stable rank mismatch");
  expect_close(spectrum.leading_spectral_gap_ratio, 4.0,
               "spectral gap mismatch");
  if (!(spectrum.normalized_spectral_entropy > 0.0 &&
        spectrum.normalized_spectral_entropy < 1.0)) {
    throw std::runtime_error("normalized entropy outside (0, 1)");
  }
  return 0;
}
