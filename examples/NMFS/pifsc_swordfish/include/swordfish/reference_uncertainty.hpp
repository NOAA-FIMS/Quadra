/** @file reference_uncertainty.hpp
 * Full correlated local delta-method uncertainty for fitted reference points.
 * Curvature uses centered differences of exact AD gradients; root-derived
 * output Jacobians use centered differences at two steps. Not bootstrap
 * uncertainty.
 */
#pragma once
#include "core/autodiff/model_gradient.hpp"
#include "local_curvature.hpp"
#include "model.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <ostream>

namespace swordfish {
/// One reference-point quantity with native units; proportions get bounded
/// intervals.
struct ReferenceEstimate {
  std::string benchmark, quantity, units;
  double value;
  bool proportion;
};
/// Flatten all six standard quantities for each benchmark into a stable 24-row
/// order. Known unfished zeros/SPR=1 and target SPR=.2 are exact constants.
/// Btgt20 raw SPR is .2/h under the supported positive Beverton-Holt
/// equilibrium, avoiding spurious root-rounding derivatives when h is fixed.
inline std::vector<ReferenceEstimate>
reference_estimates(const SwordfishModel &model, const std::vector<double> &x) {
  const auto points = model.reference_points(x);
  const auto native = model.decode(x);
  std::vector<ReferenceEstimate> result;
  const std::array<std::string, 6> names{
      "F_multiplier",      "SPR",
      "SSB_tonnes",        "recruitment_thousands",
      "dead_catch_tonnes", "annual_F_SS3_units4_ages1_10"};
  const std::array<std::string, 6> units{"multiplier", "proportion",
                                         "tonnes",     "thousands_of_fish",
                                         "tonnes",     "per_year"};
  for (const auto &point : points) {
    const auto &r = point.second;
    std::array<double, 6> values{r.scale,       r.spr,   r.ssb,
                                 r.recruitment, r.yield, r.annual_f};
    if (point.first == "unfished") {
      values[0] = 0.;
      values[1] = 1.;
      values[4] = 0.;
      values[5] = 0.;
    }
    if (point.first == "SPR20")
      values[1] = .2;
    if (point.first == "Btgt20")
      values[1] = .2 / native.first.at("SR_BH_steep");
    for (size_t i = 0; i < names.size(); ++i)
      result.push_back({point.first, names[i], units[i], values[i], i == 1});
  }
  return result;
}
/// Write a finite labeled matrix, retaining full coordinate/output
/// correlations.
inline void write_uncertainty_matrix(const std::filesystem::path &path,
                                     const Eigen::MatrixXd &matrix,
                                     const std::vector<std::string> &rows,
                                     const std::vector<std::string> &columns) {
  if (matrix.rows() != int(rows.size()) ||
      matrix.cols() != int(columns.size()) || !matrix.allFinite())
    throw std::invalid_argument("Uncertainty matrix dimensions/values");
  std::ofstream out(path);
  out << "quantity";
  for (const auto &name : columns)
    out << ',' << name;
  out << '\n' << std::setprecision(17);
  for (size_t i = 0; i < rows.size(); ++i) {
    out << rows[i];
    for (size_t j = 0; j < columns.size(); ++j)
      out << ',' << matrix(i, j);
    out << '\n';
  }
  out.flush();
  if (!out)
    throw std::runtime_error("Cannot write uncertainty matrix");
}
/// Calculate fresh curvature and propagate J V J-transpose for this fitted
/// model. The caller must have verified convergence. Every estimate uses the
/// same native parameter/deviation decoding and profiled Q as the fit. Hessian
/// and Jacobian checks at 1e-4 and 5e-5 reject numerical instability. Outputs
/// are conditional local parameter intervals; no forecast process, bootstrap or
/// model uncertainty.
inline void write_reference_point_uncertainty(
    SwordfishModel &model, const std::vector<double> &x,
    const std::filesystem::path &folder, std::ostream &progress) {
  const int n = x.size();
  const auto estimates = reference_estimates(model, x);
  const int m = estimates.size();
  std::vector<std::string> labels;
  for (const auto &r : estimates)
    labels.push_back(r.benchmark + "/" + r.quantity);
  const auto coordinates = model.parameter_names_impl();
  std::array<Eigen::MatrixXd, 2> hessian{Eigen::MatrixXd(n, n),
                                         Eigen::MatrixXd(n, n)},
      jacobian{Eigen::MatrixXd(m, n), Eigen::MatrixXd(m, n)};
  for (int pass = 0; pass < 2; ++pass) {
    const double step = pass == 0 ? 1e-4 : 5e-5;
    for (int j = 0; j < n; ++j) {
      auto plus = x, minus = x;
      plus[j] += step;
      minus[j] -= step;
      const auto gp = quadra::evaluate_gradient(model, plus),
                 gm = quadra::evaluate_gradient(model, minus);
      for (int i = 0; i < n; ++i)
        hessian[pass](i, j) =
            (gp.gradient_m[i] - gm.gradient_m[i]) / (2 * step);
      const auto a = reference_estimates(model, plus),
                 b = reference_estimates(model, minus);
      for (int i = 0; i < m; ++i)
        jacobian[pass](i, j) = (a[i].value - b[i].value) / (2 * step);
      if (j % 7 == 0 || j == n - 1)
        progress << "Reference-point uncertainty step " << step << ": column "
                 << j + 1 << '/' << n << std::endl;
    }
    const std::string suffix = pass == 0 ? "coarse" : "fine";
    write_uncertainty_matrix(folder / ("reference_hessian_" + suffix + ".csv"),
                             hessian[pass], coordinates, coordinates);
    write_uncertainty_matrix(folder / ("reference_jacobian_" + suffix + ".csv"),
                             jacobian[pass], labels, coordinates);
  }
  const auto curvature = validated_local_curvature(hessian[0], hessian[1]);
  const auto covariance = delta_covariance(jacobian[1], curvature.covariance);
  const auto coarse_j = delta_covariance(jacobian[0], curvature.covariance);
  const auto coarse_h =
      delta_covariance(jacobian[1], curvature.coarse_covariance);
  double se_change = 0.;
  for (int i = 0; i < m; ++i) {
    if (covariance(i, i) < 0 || coarse_j(i, i) < 0 || coarse_h(i, i) < 0)
      throw std::runtime_error("Negative reference-point variance");
    const double se = std::sqrt(covariance(i, i)),
                 denominator =
                     std::max(se, std::abs(estimates[i].value) * 1e-10);
    if (denominator > 0)
      se_change = std::max(
          {se_change, std::abs(std::sqrt(coarse_j(i, i)) - se) / denominator,
           std::abs(std::sqrt(coarse_h(i, i)) - se) / denominator});
  }
  if (se_change > 1e-3)
    throw std::runtime_error(
        "Reference-point SEs sensitive to numerical steps");
  write_uncertainty_matrix(folder / "reference_coordinate_covariance.csv",
                           curvature.covariance, coordinates, coordinates);
  write_uncertainty_matrix(folder / "management_reference_point_covariance.csv",
                           covariance, labels, labels);
  std::ofstream out(folder / "management_reference_point_uncertainty.csv");
  out << "benchmark,quantity,units,estimate,standard_error,CV,lower95,upper95,"
         "interval_method\n"
      << std::setprecision(17);
  for (int i = 0; i < m; ++i) {
    const auto &r = estimates[i];
    const double se = std::sqrt(covariance(i, i));
    double lower = r.value, upper = r.value;
    std::string method = "constant";
    if (se > 0) {
      if (r.proportion) {
        if (!(r.value > 0 && r.value < 1))
          throw std::runtime_error("Invalid SPR proportion interval");
        const double center = std::log(r.value / (1 - r.value)),
                     width = 1.95996398454005 * se / (r.value * (1 - r.value));
        lower = 1 / (1 + std::exp(-(center - width)));
        upper = 1 / (1 + std::exp(-(center + width)));
        method = "logit_Wald_95";
      } else {
        const auto interval =
            positive_delta_interval(r.value, covariance(i, i));
        lower = interval.lower;
        upper = interval.upper;
        method = "log_Wald_95";
      }
    }
    out << r.benchmark << ',' << r.quantity << ',' << r.units << ',' << r.value
        << ',' << se << ',';
    if (r.value != 0)
      out << se / std::abs(r.value);
    out << ',' << lower << ',' << upper << ',' << method << '\n';
  }
  out.flush();
  if (!out)
    throw std::runtime_error("Cannot write reference-point uncertainty");
  std::ofstream status(folder / "reference_uncertainty_summary.json");
  status
      << std::setprecision(17)
      << "{\"passed\":true,\"method\":\"full correlated local delta method; "
         "joint penalized likelihood\",\"hessian\":\"central differences of "
         "exact AD gradients\",\"jacobian\":\"central differences through "
         "model reference-point "
         "solves\",\"steps\":[0.0001,0.00005],\"coordinates\":"
      << n << ",\"quantities\":" << m
      << ",\"relative_hessian_asymmetry\":" << curvature.asymmetry
      << ",\"relative_hessian_step_difference\":" << curvature.step_difference
      << ",\"minimum_eigenvalue\":" << curvature.minimum_eigenvalue
      << ",\"maximum_eigenvalue\":" << curvature.maximum_eigenvalue
      << ",\"maximum_parameter_SE_step_change\":"
      << curvature.maximum_se_step_change
      << ",\"maximum_reference_SE_step_change\":" << se_change
      << ",\"regularization_used\":false,\"saved_covariance_used\":false,"
         "\"confidence_level\":0.95,\"process_uncertainty_included\":false}\n";
  status.flush();
  if (!status)
    throw std::runtime_error("Cannot write reference uncertainty status");
  progress << "Reference-point uncertainty validated; intervals written."
           << std::endl;
}
} // namespace swordfish
