/**
 * @file swordfish_projection_uncertainty.cpp
 * @details
 * Local full-covariance propagation for five scenarios, ten years and five
 * quantities. Every perturbation recalibrates the scenario rule and fishing
 * mixture. Positive pointwise intervals omit future recruitment process
 * variability and are not assessment bootstrap intervals.
 */
#include "swordfish/delta_method.hpp"
#include "swordfish/model.hpp"
#include "swordfish/reference_io.hpp"
#include "swordfish/scenarios.hpp"
#include "swordfish_input.hpp"
#include <filesystem>
#include <iomanip>
#include <iostream>
using namespace swordfish;
/// One labeled scenario-year-quantity point estimate, with a unique
/// covariance/Jacobian output identifier.
struct Output {
  int scenario, year;
  std::string name, quantity, label;
  double estimate;
};
/// Flatten all scenario trajectories into a stable 250-output order shared by
/// intervals, Jacobian rows and covariance labels.
std::vector<Output> flatten(const std::vector<ScenarioProjection> &scenarios) {
  std::vector<Output> out;
  for (const auto &s : scenarios)
    for (const auto &r : s.years) {
      for (const auto &item : std::vector<std::pair<std::string, double>>{
               {"spawning_biomass_tonnes", r.ssb},
               {"recruitment_thousands", r.recruitment},
               {"biomass_january_tonnes", r.biomass},
               {"dead_catch_tonnes", r.yield},
               {"annual_F_SS3_units4_ages1_10", r.annual_f}})
        out.push_back({s.scenario.definition.id, r.year,
                       s.scenario.definition.name, item.first,
                       "S" + std::to_string(s.scenario.definition.id) + "_Y" +
                           std::to_string(r.year) + "_" + item.first,
                       item.second});
    }
  return out;
}
/// Run swordfish projection uncertainty; return zero on success and nonzero on
/// a failed check or invalid input.
int main() {
  try {
    auto data = imported_swordfish_data();
    auto native =
        read_fitted_parameters(data, "build/quadra_fitted_parameters.csv",
                               "build/biology_fitted_parameters.txt",
                               "build/population_fitted_deviations.txt");
    SwordfishModel model(data, native.first, native.second);
    const auto x = model.initial_vector();
    const auto definitions =
        read_projection_scenarios("data/projection_scenarios.csv");
    const int n = x.size();
    // Covariance must belong to exactly the current coordinate vector, in the
    // same order.
    std::ifstream coordinates("build/uncertainty_coordinates.csv");
    std::string line;
    if (!csv_line(coordinates, line) ||
        line != "parameter,coordinate,native,lower,upper")
      throw std::runtime_error("Missing covariance coordinates");
    for (int i = 0; i < n; ++i) {
      if (!csv_line(coordinates, line))
        throw std::runtime_error("Incomplete covariance coordinates");
      auto c = csv_fields(line);
      if (c.size() != 5 || c[0] != model.free[i].name ||
          !std::isfinite(csv_number(c[1])) ||
          std::abs(csv_number(c[1]) - x[i]) > 1e-12 ||
          csv_number(c[3]) != model.free[i].lo ||
          csv_number(c[4]) != model.free[i].hi)
        throw std::runtime_error(
            "Covariance is stale or coordinate mapping differs");
    }
    if (csv_line(coordinates, line))
      throw std::runtime_error("Extra covariance coordinates");
    std::ifstream in("build/coordinate_covariance.csv");
    if (!csv_line(in, line))
      throw std::runtime_error("Missing covariance");
    auto header = csv_fields(line);
    if (header.size() != x.size() + 1 || header[0] != "parameter")
      throw std::runtime_error("Covariance header");
    Eigen::MatrixXd covariance(n, n);
    for (int i = 0; i < n; ++i) {
      if (header[i + 1] != model.free[i].name || !csv_line(in, line))
        throw std::runtime_error("Covariance column labels");
      auto c = csv_fields(line);
      if (c.size() != x.size() + 1 || c[0] != model.free[i].name)
        throw std::runtime_error("Covariance row labels");
      for (int j = 0; j < n; ++j)
        covariance(i, j) = csv_number(c[j + 1]);
    }
    if (csv_line(in, line))
      throw std::runtime_error("Extra covariance rows");
    validate_parameter_covariance(covariance);
    auto evaluate = [&](const std::vector<double> &point) {
      const auto decoded = model.decode(point);
      return flatten(evaluate_projection_scenarios(
          data, decoded.first, decoded.second, definitions));
    };
    const auto center = evaluate(x);
    const int m = center.size();
    if (m != 250)
      throw std::runtime_error("Expected 250 projected outputs");
    auto values = [&](const std::vector<double> &point) {
      const auto result = evaluate(point);
      if (result.size() != center.size())
        throw std::runtime_error("Projection output size changed");
      Eigen::VectorXd out(m);
      for (int i = 0; i < m; ++i) {
        if (result[i].label != center[i].label || !(result[i].estimate > 0) ||
            !std::isfinite(result[i].estimate))
          throw std::runtime_error("Invalid projected value or order");
        out[i] = result[i].estimate;
      }
      return out;
    };
    const std::filesystem::path folder = "build/projections";
    std::filesystem::create_directories(folder);
    // Rows are scenario/year/quantity outputs; columns are the 84 independent
    // logit coordinates. Central differences include the whole composed
    // mapping: historical state -> reference pattern/target -> projected
    // population.
    Eigen::MatrixXd coarse(m, n), fine(m, n);
    for (int pass = 0; pass < 2; ++pass) {
      const double step = pass == 0 ? 1e-4 : 5e-5;
      auto &jac = pass == 0 ? coarse : fine;
      for (int j = 0; j < n; ++j) {
        auto plus = x, minus = x;
        plus[j] += step;
        minus[j] -= step;
        jac.col(j) = (values(plus) - values(minus)) / (2 * step);
        if (j % 14 == 0)
          std::cout << "Projection derivatives, step " << step
                    << ", coordinate " << j << '/' << n << std::endl;
      }
    }
    // Never substitute a diagonal covariance here. Historical abundance and F
    // estimates are correlated, and target recalibration can compensate
    // effects. The 250-output covariance can be singular because there are only
    // 84 inputs.
    const Eigen::MatrixXd propagated = delta_covariance(fine, covariance),
                          coarse_cov = delta_covariance(coarse, covariance);
    double max_se_change = 0, max_gradient_change = 0;
    for (int i = 0; i < m; ++i) {
      if (!(propagated(i, i) > 0))
        throw std::runtime_error("Nonpositive projected variance");
      max_se_change = std::max(
          max_se_change,
          std::abs(std::sqrt(coarse_cov(i, i) / propagated(i, i)) - 1));
      max_gradient_change =
          std::max(max_gradient_change,
                   (coarse.row(i) - fine.row(i)).norm() /
                       std::max(fine.row(i).norm(), 1e-8 * center[i].estimate));
    }
    // Independent directions test the entire composed target and trajectory
    // calculation.
    double max_direction_error = 0;
    for (int direction = 0; direction < 4; ++direction) {
      Eigen::VectorXd v(n);
      for (int j = 0; j < n; ++j)
        v[j] = std::sin((j + 1.) * (direction + 1.37));
      v.normalize();
      auto plus = x, minus = x;
      for (int j = 0; j < n; ++j) {
        plus[j] += 1e-4 * v[j];
        minus[j] -= 1e-4 * v[j];
      }
      const Eigen::VectorXd numeric = (values(plus) - values(minus)) / 2e-4,
                            analytic = fine * v;
      for (int i = 0; i < m; ++i)
        max_direction_error = std::max(
            max_direction_error,
            std::abs(numeric[i] - analytic[i]) /
                std::max(std::abs(analytic[i]), 1e-3 * center[i].estimate));
    }
    if (max_se_change > 1e-4 || max_gradient_change > 1e-3 ||
        max_direction_error > 1e-3)
      throw std::runtime_error(
          "Projection derivatives are sensitive to step or direction");
    std::ofstream out(folder / "projection_uncertainty.csv"),
        jacobian(folder / "projection_jacobian.csv"),
        cov_out(folder / "projection_output_covariance.csv"),
        summary(folder / "projection_uncertainty_numerics.json");
    out << "scenario,name,year,quantity,estimate,standard_error,coefficient_of_"
           "variation,lower95,upper95\n"
        << std::setprecision(17);
    for (int i = 0; i < m; ++i) {
      const auto &c = center[i];
      const auto interval =
          positive_delta_interval(c.estimate, propagated(i, i));
      out << c.scenario << ',' << c.name << ',' << c.year << ',' << c.quantity
          << ',' << c.estimate << ',' << interval.se << ','
          << interval.se / c.estimate << ',' << interval.lower << ','
          << interval.upper << '\n';
    }
    jacobian << "output";
    for (const auto &p : model.free)
      jacobian << ',' << p.name;
    jacobian << '\n' << std::setprecision(17);
    cov_out << "output";
    for (const auto &c : center)
      cov_out << ',' << c.label;
    cov_out << '\n' << std::setprecision(17);
    for (int i = 0; i < m; ++i) {
      jacobian << center[i].label;
      for (int j = 0; j < n; ++j)
        jacobian << ',' << fine(i, j);
      jacobian << '\n';
      cov_out << center[i].label;
      for (int j = 0; j < m; ++j)
        cov_out << ',' << propagated(i, j);
      cov_out << '\n';
    }
    summary
        << std::setprecision(17) << "{\"passed\":true,\"parameters\":" << n
        << ",\"outputs\":" << m
        << ",\"steps\":[0.0001,0.00005],\"maximum_relative_SE_step_change\":"
        << max_se_change
        << ",\"maximum_relative_gradient_step_change\":" << max_gradient_change
        << ",\"directional_checks\":4,\"maximum_scaled_directional_error\":"
        << max_direction_error
        << ",\"interval\":\"pointwise 95% log-scale Wald; native SE via delta "
           "method\",\"target_rule\":\"recalibrate scenario targets and "
           "fishing pattern for every parameter "
           "perturbation\",\"future_recruitment_process_variation\":false}\n";
    if (!out || !jacobian || !cov_out || !summary)
      throw std::runtime_error("Projection uncertainty output failed");
    std::cout << "Wrote 250 projection intervals; maximum SE step change "
              << max_se_change << std::endl;
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
