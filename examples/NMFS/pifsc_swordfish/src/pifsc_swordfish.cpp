/** @file pifsc_swordfish.cpp
 * Standalone C++ SS3 conversion, JSON intake, evaluation, fitting and
 * reporting. The private stdin transport remains available for existing
 * analysis scripts.
 */
#include "core/autodiff/model_gradient.hpp"
#include "core/optimizer.hpp"
#include "swordfish/csv.hpp"
#include "swordfish/model.hpp"
#include "swordfish/native_cli.hpp"
#include "swordfish/native_io.hpp"
#include "swordfish/optimizer_scaling.hpp"
#include "swordfish/reference_uncertainty.hpp"
#include "swordfish/report.hpp"
#include "swordfish/ss3_native.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
DECLARE_ADGRAPH();
namespace {
void run_model(
    std::pair<swordfish::InputData, swordfish::Deviations<double>> native,
    const swordfish::NativeOptions &options, bool &created) {
  const auto &output = options.output;
  const auto &mode = options.mode;
  const auto &metric_path = options.metric;
  const auto iterations = options.iterations;
  const auto tolerance = options.tolerance;
  const auto uncertainty = options.uncertainty;
  auto parameters = swordfish::initial_values<double>(native.first);
  swordfish::SwordfishModel model(
      std::move(native.first), std::move(parameters), std::move(native.second));
  auto x = model.initial_vector();
  const auto initial = quadra::evaluate_gradient(model, x);
  if (!std::isfinite(initial.objective_value_m) ||
      !std::isfinite(initial.gradient_norm_m))
    throw std::invalid_argument("Nonfinite initial objective/gradient");
  if (!output.parent_path().empty())
    std::filesystem::create_directories(output.parent_path());
  if (!std::filesystem::create_directory(output))
    throw std::invalid_argument("Output directory already exists");
  created = true;
  std::cout << std::setprecision(15) << "Native model: " << x.size()
            << " independent coordinates; initial objective "
            << initial.objective_value_m << std::endl;
  bool converged = false;
  int completed = 0;
  if (mode == "fit") {
    quadra::OptResult fit;
    if (metric_path.empty()) {
      auto initial_parameters = model.initial_parameters();
      fit = quadra::optimize_lbfgs(model, initial_parameters,
                                   quadra::default_laplace_options(),
                                   iterations, tolerance);
      x = fit.par;
    } else {
      const auto covariance = swordfish::read_optimizer_metric(
          metric_path, model.parameter_names_impl());
      Eigen::LLT<Eigen::MatrixXd> decomposition(covariance);
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(covariance);
      if (decomposition.info() != Eigen::Success ||
          eigen.info() != Eigen::Success)
        throw std::runtime_error("Optimizer metric factorization failed");
      swordfish::OptimizerScaling<swordfish::SwordfishModel> scaled(
          model, x, decomposition.matrixL());
      quadra::ParameterVector initial_parameters;
      for (const auto &name : model.parameter_names_impl())
        initial_parameters.add(quadra::Parameter(name, 0.));
      const double scaled_tolerance = std::min(
          1e-7, tolerance * std::sqrt(eigen.eigenvalues().minCoeff()) * .1);
      fit = quadra::optimize_lbfgs(scaled, initial_parameters,
                                   quadra::default_laplace_options(),
                                   iterations, scaled_tolerance);
      x = scaled.expand(fit.par);
    }
    std::ofstream trace(output / "quadra_fit_trace.csv");
    trace << "evaluation,objective,gradient_norm\n" << std::setprecision(17);
    for (const auto &row : fit.trajectory)
      trace << row.evaluation << ',' << row.objective << ','
            << row.gradient_norm << '\n';
    trace.flush();
    if (!trace)
      throw std::runtime_error("Cannot save native fit trace");
    if (!fit.converged)
      throw std::runtime_error(
          "Native fit did not converge; fitted outputs were not written");
    completed = fit.iterations;
    converged = true;
  }
  const auto checked = quadra::evaluate_gradient(model, x);
  if (!std::isfinite(checked.objective_value_m) ||
      !std::isfinite(checked.gradient_norm_m) ||
      (mode == "fit" && checked.gradient_norm_m > tolerance))
    throw std::runtime_error("Native final gradient check failed");
  auto fitted = model.decode(x);
  const auto likelihood = model.likelihood(x);
  for (const auto &row : likelihood.indices)
    fitted.first.at("LnQ_fleet_" + std::to_string(row.fleet)) = row.log_q;
  swordfish::write_fitted_outputs(model.data, fitted.first, fitted.second,
                                  output);
  std::ofstream values(output / "quadra_fitted_parameters.csv");
  values << "parameter,value\n" << std::setprecision(17);
  for (const auto &row : fitted.first)
    values << row.first << ',' << row.second << '\n';
  for (const auto &row : fitted.second)
    values << "RecruitmentDev_" << row.first << ',' << row.second << '\n';
  std::ofstream gradient(output / "objective_gradient.csv");
  gradient << "parameter,unconstrained_value,gradient\n"
           << std::setprecision(17);
  for (size_t i = 0; i < x.size(); ++i)
    gradient << model.free[i].name << ',' << x[i] << ','
             << checked.gradient_m[i] << '\n';
  if (uncertainty)
    swordfish::write_reference_point_uncertainty(model, x, output, std::cout);
  std::ofstream summary(output / "quadra_fit_summary.json");
  summary << std::setprecision(17) << "{\"mode\":\"" << mode
          << "\",\"initialization\":\"native normalized JSON\",\"objective\":"
          << checked.objective_value_m
          << ",\"initial_objective\":" << initial.objective_value_m
          << ",\"gradient_norm\":" << checked.gradient_norm_m
          << ",\"converged\":" << (converged ? "true" : "false")
          << ",\"free_parameters\":" << x.size()
          << ",\"optimizer_iterations\":" << completed
          << ",\"optimizer_metric_used\":"
          << (metric_path.empty() ? "false" : "true")
          << ",\"reference_points_calculated\":true,\"reference_uncertainty_"
             "calculated\":"
          << (uncertainty ? "true" : "false") << "}\n";
  values.flush();
  gradient.flush();
  summary.flush();
  if (!values || !gradient || !summary)
    throw std::runtime_error("Native output write failure");
  std::cout << "Objective: " << checked.objective_value_m
            << "; gradient norm: " << checked.gradient_norm_m
            << "; reference points written." << std::endl;
}
} // namespace

/// Convert SS3 inputs or run the PIFSC swordfish model directly from JSON.
int main(int argc, char **argv) {
  swordfish::NativeOptions options;
  bool created = false;
  std::string source, metric_hash, binary_hash;
  int result = 0;
  std::string error_message;
  try {
    options = swordfish::native_options(argc, argv);
    if (options.help) {
      swordfish::print_native_help(std::cout, options.convert);
      return 0;
    }
    if (options.convert) {
      swordfish::convert_ss3(options.data, options.control, options.output,
                             options.parameters);
      return 0;
    }
    swordfish::require_new_path(options.output);
    auto native = [&]() {
      if (options.legacy)
        return swordfish::read_native_input(std::cin);
      source = swordfish::read_file(options.input);
      return swordfish::native_from_json(swordfish::parse_json(source));
    }();
    if (!options.legacy) {
      binary_hash =
          swordfish::sha256(swordfish::read_file(swordfish::executable_path()));
      if (!options.metric.empty())
        metric_hash = swordfish::sha256(swordfish::read_file(options.metric));
      std::cout << "Native input: "
                << swordfish::fs::absolute(options.input).string() << std::endl;
    }
    run_model(std::move(native), options, created);
  } catch (const swordfish::CLIError &error) {
    std::cerr << "Native command error: " << error.what() << '\n';
    return 2;
  } catch (const std::exception &error) {
    result = 1;
    error_message = error.what();
    std::cerr << "Native model error: " << error_message << '\n';
  }
  if (created && !options.legacy) {
    bool metric_unchanged = true;
    try {
      if (!options.metric.empty())
        metric_unchanged =
            swordfish::sha256(swordfish::read_file(options.metric)) ==
            metric_hash;
    } catch (const std::exception &) {
      metric_unchanged = false;
    }
    if (!metric_unchanged) {
      result = 1;
      error_message = "Optimizer metric changed during run";
      std::cerr << "Native model error: " << error_message << '\n';
    }
    try {
      swordfish::exclusive_file(options.output / "model_input.json", source);
      swordfish::Json status = {
          {"status", result == 0 ? "completed" : "failed"},
          {"mode", options.mode},
          {"input_sha256", swordfish::sha256(source)},
          {"binary_sha256", binary_hash},
          {"ss3_files_used", false},
          {"uncertainty_requested", options.uncertainty},
          {"optimizer_metric_sha256", metric_hash.empty()
                                          ? swordfish::Json(nullptr)
                                          : swordfish::Json(metric_hash)},
          {"optimizer_metric_unchanged", metric_unchanged},
          {"returncode", result}};
      if (!error_message.empty())
        status["error"] = error_message;
      swordfish::exclusive_file(options.output / "run_status.json",
                                status.dump(2) + '\n');
    } catch (const std::exception &error) {
      std::cerr << "Native provenance error: " << error.what() << '\n';
      result = 1;
    }
  }
  return result;
}
