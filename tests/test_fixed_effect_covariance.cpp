#include "../core/inference/fixed_effect_covariance.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

struct QuadraticObjective {
  Eigen::MatrixXd H;

  double operator()(const std::vector<double> &x) const {
    Eigen::VectorXd v(static_cast<Eigen::Index>(x.size()));
    for (Eigen::Index i = 0; i < v.size(); ++i) {
      v[i] = x[static_cast<size_t>(i)];
    }

    return 0.5 * v.transpose() * H * v;
  }
};

struct QuadraticGradient {
  Eigen::MatrixXd H;
  int calls = 0;

  std::vector<double> operator()(const std::vector<double> &x) {
    ++calls;
    const Eigen::Map<const Eigen::VectorXd> v(
        x.data(), static_cast<Eigen::Index>(x.size()));
    const Eigen::VectorXd gradient = H * v;
    return std::vector<double>(gradient.data(),
                               gradient.data() + gradient.size());
  }
};

int main() {
  QuadraticObjective objective;
  objective.H.resize(2, 2);
  objective.H << 4.0, 1.0, 1.0, 3.0;

  const std::vector<double> theta_hat{0.2, -0.3};

  auto result =
      quadra::estimate_fixed_effect_covariance(objective, theta_hat, 1.0e-4);

  if (!result.success_m) {
    std::cerr << "FAIL: covariance estimation failed: " << result.message_m
              << "\n";
    return 1;
  }

  const Eigen::MatrixXd expected_covariance = objective.H.inverse();

  const auto direct_result =
      quadra::estimate_fixed_effect_covariance_from_hessian(objective.H);
  if (!direct_result.success_m ||
      (direct_result.covariance_m - expected_covariance)
              .cwiseAbs()
              .maxCoeff() > 1.0e-12) {
    std::cerr << "FAIL: direct-Hessian covariance path is inaccurate\n";
    return 1;
  }

  const double hessian_error =
      (result.hessian_m - objective.H).cwiseAbs().maxCoeff();

  const double covariance_error =
      (result.covariance_m - expected_covariance).cwiseAbs().maxCoeff();

  if (hessian_error > 1.0e-6) {
    std::cerr << "FAIL: Hessian error too large: " << hessian_error << "\n";
    return 1;
  }

  if (covariance_error > 1.0e-6) {
    std::cerr << "FAIL: covariance error too large: " << covariance_error
              << "\n";
    return 1;
  }

  if (std::abs(result.correlation_m(0, 0) - 1.0) > 1.0e-12 ||
      std::abs(result.correlation_m(1, 1) - 1.0) > 1.0e-12) {
    std::cerr << "FAIL: correlation diagonal is not 1.\n";
    return 1;
  }

  QuadraticGradient gradient{objective.H};
  const auto gradient_result =
      quadra::estimate_fixed_effect_covariance_from_gradient(gradient,
                                                             theta_hat, 1.0e-4);
  if (!gradient_result.success_m || gradient.calls != 4) {
    std::cerr << "FAIL: gradient covariance path did not use exactly 2p "
                 "evaluations\n";
    return 1;
  }
  if ((gradient_result.hessian_m - objective.H).cwiseAbs().maxCoeff() >
          1.0e-10 ||
      (gradient_result.covariance_m - expected_covariance)
              .cwiseAbs()
              .maxCoeff() > 1.0e-10) {
    std::cerr << "FAIL: gradient covariance path is inaccurate\n";
    return 1;
  }

  std::cout << "PASS: fixed-effect covariance finite-difference test\n";
  std::cout << "  max Hessian error: " << hessian_error << "\n";
  std::cout << "  max covariance error: " << covariance_error << "\n";
  std::cout << "  gradient evaluations: " << gradient.calls << "\n";

  return 0;
}
