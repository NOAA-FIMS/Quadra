/**
 * @file delta_method.hpp
 * @brief Full correlated local parameter uncertainty.
 *
 * Input covariance is for independent optimizer coordinates. Propagate J V
 * J-transpose with all correlations retained; output covariance can be singular
 * when outputs outnumber parameters and must not be inverted as if positive
 * definite. Log-scale Wald endpoints are pointwise local approximations, not
 * simultaneous bands or future-realization prediction intervals.
 */
#pragma once
#include <Eigen/Dense>
#include <cmath>
#include <stdexcept>
namespace swordfish {
/// Require a nonempty finite, symmetric, positive-definite parameter covariance
/// before propagation; do not regularize silently.
/// @param covariance Square, finite, symmetric positive-definite covariance of
/// independent fitting coordinates.
inline void validate_parameter_covariance(const Eigen::MatrixXd &covariance) {
  if (covariance.rows() == 0 || covariance.rows() != covariance.cols() ||
      !covariance.allFinite())
    throw std::invalid_argument(
        "Invalid parameter covariance dimensions or values");
  if ((covariance - covariance.transpose()).norm() > 1e-10 * covariance.norm())
    throw std::invalid_argument("Asymmetric parameter covariance");
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(covariance);
  if (solver.info() != Eigen::Success || solver.eigenvalues().minCoeff() <= 0)
    throw std::invalid_argument(
        "Parameter covariance must be positive definite");
}
/// Propagate a compatible finite output Jacobian through the full parameter
/// covariance as J V J-transpose, preserving correlations.
/// @param jacobian Output-by-coordinate derivative matrix, with columns in the
/// covariance coordinate order.
/// @param covariance Square, finite, symmetric positive-definite covariance of
/// independent fitting coordinates.
/// @return Native output covariance J*V*J-transpose, which may be rank
/// deficient.
/// @details Columns of jacobian must correspond exactly to covariance rows.
/// Keeping only marginal parameter variances would omit compensating/correlated
/// effects and change the reported uncertainties.
inline Eigen::MatrixXd delta_covariance(const Eigen::MatrixXd &jacobian,
                                        const Eigen::MatrixXd &covariance) {
  // Validate before multiplication so Eigen assertions are not the public
  // error path for an invalid matrix. Never silently repair an indefinite fit.
  validate_parameter_covariance(covariance);
  if (jacobian.cols() != covariance.rows() || !jacobian.allFinite())
    throw std::invalid_argument("Invalid derived Jacobian");
  return jacobian * covariance * jacobian.transpose();
}
/// Native standard error and positive lower/upper pointwise log-scale Wald
/// endpoints.
/// @details se is native-scale delta-method uncertainty. lower/upper are
/// positive log-scale Wald endpoints, so their arithmetic midpoint need not
/// equal the estimate. No process uncertainty is added.
struct PositiveInterval {
  double se, lower, upper;
};
/// Construct a pointwise positive 95% interval from a positive estimate and
/// nonnegative native variance; reject nonfinite results.
/// @param estimate Strictly positive native-scale output estimate.
/// @param variance Nonnegative finite native-scale marginal variance, not a
/// standard deviation.
/// @return Native SE and positive pointwise 95% log-scale Wald bounds.
/// @details Bounds are estimate*exp(+-1.95996398454005*SE/estimate). They are
/// not bootstrap quantiles, simultaneous trajectory bands or intervals for
/// realized future recruitment.
inline PositiveInterval positive_delta_interval(double estimate,
                                                double variance) {
  if (!(estimate > 0) || !std::isfinite(estimate) || !(variance >= 0) ||
      !std::isfinite(variance))
    throw std::invalid_argument("Positive delta interval inputs");
  const double se = std::sqrt(variance),
               width = 1.95996398454005 * se / estimate;
  const double lower = estimate * std::exp(-width),
               upper = estimate * std::exp(width);
  if (!(lower > 0) || !std::isfinite(upper))
    throw std::runtime_error("Nonfinite local interval");
  return {se, lower, upper};
}
} // namespace swordfish
