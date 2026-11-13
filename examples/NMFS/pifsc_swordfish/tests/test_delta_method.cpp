#include "swordfish/delta_method.hpp"
#include <iostream>
/// Test helper for test delta method: near. See the typed signature and
/// assertions in the linked test source.
void near(double a, double b) {
  if (std::abs(a - b) > 1e-10 * std::max(1., std::abs(b)))
    throw std::runtime_error("Delta-method analytic check failed");
}
/// Run test delta method; return zero on success and nonzero on a failed check
/// or invalid input.
int main() {
  try {
    Eigen::Matrix2d sigma;
    sigma << 4, -1, -1, 9;
    Eigen::Matrix2d jac;
    jac << 1, 2, -3, 1;
    swordfish::validate_parameter_covariance(sigma);
    const auto result = swordfish::delta_covariance(jac, sigma);
    near(result(0, 0), 36);
    near(result(1, 1), 51);
    near(result(0, 1), 11);
    near(result(1, 0), 11);
    // A linear change of parameter coordinates must leave all propagated
    // covariances unchanged.
    Eigen::Matrix2d scale = Eigen::Matrix2d::Zero();
    scale(0, 0) = 2;
    scale(1, 1) = .5;
    const auto transformed = swordfish::delta_covariance(
        jac * scale, scale.inverse() * sigma * scale.inverse().transpose());
    if ((transformed - result).norm() > 1e-10)
      throw std::runtime_error("Coordinate invariance failure");
    const auto interval = swordfish::positive_delta_interval(10, 36);
    near(interval.se, 6);
    near(std::sqrt(interval.lower * interval.upper), 10);
    const auto zero = swordfish::positive_delta_interval(10, 0);
    near(zero.lower, 10);
    near(zero.upper, 10);
    bool rejected = false;
    try {
      Eigen::Matrix2d invalid;
      invalid << 1, 2, 2, 1;
      swordfish::validate_parameter_covariance(invalid);
    } catch (const std::invalid_argument &) {
      rejected = true;
    }
    if (!rejected)
      throw std::runtime_error("Indefinite covariance accepted");
    std::cout << "Correlated delta covariance, coordinate invariance, and "
                 "interval checks pass\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
