/**
 * @file swordfish_control_start.cpp
 * @brief Fit from control initial values without reading Report.sso or fitted
 * covariance.
 *
 * Native parameters come from the original control factory; all 57 historical
 * log recruitment deviations start at zero. Activation phases use control phase
 * flags, but all stages evaluate the same final weighted joint objective. This
 * is a Quadra convergence experiment, not an exact replay of SS3's optimizer.
 */
#include "core/autodiff/model_gradient.hpp"
#include "core/optimizer.hpp"
#include "swordfish/control_start.hpp"
#include "swordfish/csv.hpp"
#include "swordfish/report.hpp"
#include "swordfish_control_input.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
DECLARE_ADGRAPH();

/// Optimize in a constant metric learned exclusively from the current cold
/// trajectory. No statistical Jacobian is added: this is a linear MLE
/// coordinate change.
struct ColdMetricModel : quadra::QuadraModel<ColdMetricModel> {
  swordfish::SwordfishModel &model;
  std::vector<double> center;
  Eigen::MatrixXd factor;
  /// Borrow the objective and copy the cold checkpoint and positive
  /// optimization metric.
  ColdMetricModel(swordfish::SwordfishModel &m, std::vector<double> point,
                  Eigen::MatrixXd matrix)
      : model(m), center(std::move(point)), factor(std::move(matrix)) {}
  /// Preserve the canonical coordinate label ordering for the public optimizer.
  std::vector<std::string> parameter_names_impl() const {
    return model.parameter_names_impl();
  }
  /// Map scaled displacements back to the original bounded-logit coordinates.
  template <class T> std::vector<T> expand(const std::vector<T> &z) const {
    if (z.size() != center.size())
      throw std::invalid_argument("Cold metric dimension");
    std::vector<T> x;
    for (size_t i = 0; i < z.size(); ++i) {
      T value = T(center[i]);
      for (size_t j = 0; j < z.size(); ++j)
        value += T(factor(i, j)) * z[j];
      x.push_back(value);
    }
    return x;
  }
  /// Evaluate the unchanged full joint objective at the mapped coordinates.
  template <class T> T operator()(const std::vector<T> &z) const {
    return model(expand(z));
  }
  /// Provide the same scalar target through Quadra's model interface.
  template <class T>
  T evaluate_impl(const std::vector<T> &z, quadra::ModelReportContext &) const {
    return (*this)(z);
  }
};

/// Write native values and all deviations; this path never changes the
/// canonical fit.
void save_parameters(const std::filesystem::path &path,
                     const swordfish::ParameterValues<double> &p,
                     const swordfish::Deviations<double> &dev) {
  std::ofstream out(path);
  out << "parameter,value\n" << std::setprecision(17);
  for (const auto &row : p)
    out << row.first << ',' << row.second << '\n';
  for (const auto &row : dev)
    out << "RecruitmentDev_" << row.first << ',' << row.second << '\n';
  if (!out)
    throw std::runtime_error("Cannot save control-start parameters");
}

/// Run strict controls or one explicit fixed-growth matching sensitivity.
/// A new output directory is mandatory so failed attempts remain reviewable.
int main(int argc, char **argv) {
  try {
    if (argc < 3 || argc > 5)
      throw std::invalid_argument(
          "Usage: swordfish_control_start strict|matched-fixed NEW_OUTPUT "
          "[refinement_iterations_per_pass] [COLD_CHECKPOINT_DIRECTORY]");
    const std::string mode = argv[1];
    if (mode != "strict" && mode != "matched-fixed")
      throw std::invalid_argument("Control-start mode");
    const std::filesystem::path folder = argv[2];
    if (std::filesystem::exists(folder))
      throw std::runtime_error("Control-start output already exists");
    std::filesystem::create_directories(folder);
    auto data = imported_swordfish_control_data();
    auto seed = swordfish::control_start_values(data, mode == "matched-fixed");
    auto p = seed.first;
    auto dev = seed.second;
    save_parameters(folder / "initial_parameters.csv", p, dev);
    swordfish::SwordfishModel model(data, p, dev);
    auto x = model.initial_vector();
    double initial = 0;
    bool initial_valid = true;
    try {
      initial = model(x);
    } catch (const std::invalid_argument &error) {
      if (std::string(error.what()) != "Empty composition")
        throw;
      initial_valid = false;
    }
    std::ofstream initial_status(folder / "initial_status.json");
    initial_status << "{\"full_objective_evaluable\":"
                   << (initial_valid ? "true" : "false") << ",\"error\":\""
                   << (initial_valid ? "" : "Empty composition") << "\"}\n";
    initial_status.close();
    std::cout << std::setprecision(15)
              << "Literal control start evaluable: " << initial_valid
              << "; mode " << mode << std::endl;
    std::ofstream trace(folder / "fit_trace.csv"),
        stages(folder / "stages.csv");
    trace << "phase,evaluation,objective,active_gradient_norm\n"
          << std::setprecision(17);
    stages << "phase,active_coordinates,initial_objective,final_objective,"
              "active_gradient_norm,converged,optimizer_iterations\n"
           << std::setprecision(17);
    const int limit = argc >= 4 ? swordfish::csv_integer(argv[3]) : 1500;
    if (limit < 1)
      throw std::invalid_argument("Positive iteration limit required");
    const std::string resume = argc == 5 ? argv[4] : "";
    if (!resume.empty()) {
      // Resume only a native checkpoint recorded by this cold experiment.
      // Seed and fixed-value checks reject incompatible configurations.
      const auto source = std::filesystem::path(resume);
      auto checkpoint_path = source / "phase_4_parameters.csv";
      for (int pass = 0; pass < 4; ++pass)
        if (std::filesystem::exists(
                source /
                ("refinement_" + std::to_string(pass) + "_parameters.csv")))
          checkpoint_path = source / ("refinement_" + std::to_string(pass) +
                                      "_parameters.csv");
      std::ifstream initial_file(source / "initial_parameters.csv"),
          checkpoint_file(checkpoint_path);
      auto read = [&](std::ifstream &in) {
        auto values = p;
        auto deviations = dev;
        std::string line;
        std::set<std::string> seen;
        if (!swordfish::csv_line(in, line) || line != "parameter,value")
          throw std::runtime_error("Cold checkpoint header");
        while (swordfish::csv_line(in, line)) {
          auto cells = swordfish::csv_fields(line);
          if (cells.size() != 2 || !seen.insert(cells[0]).second)
            throw std::runtime_error("Cold checkpoint record");
          double value = swordfish::csv_number(cells[1]);
          if (cells[0].rfind("RecruitmentDev_", 0) == 0)
            deviations.at(swordfish::csv_integer(cells[0].substr(15))) = value;
          else
            values.at(cells[0]) = value;
        }
        if (seen.size() != p.size() + dev.size())
          throw std::runtime_error("Incomplete cold checkpoint");
        return std::make_pair(values, deviations);
      };
      const auto original = read(initial_file);
      if (original.first != p || original.second != dev)
        throw std::runtime_error(
            "Resume source does not have the declared literal control seed");
      const auto checkpoint = read(checkpoint_file);
      for (const auto &spec : data.parameters)
        if (spec.columns[6] <= 0 &&
            checkpoint.first.at(spec.name) != p.at(spec.name))
          throw std::runtime_error("Resume changed a fixed control value");
      swordfish::SwordfishModel resumed(data, checkpoint.first,
                                        checkpoint.second);
      x = resumed.initial_vector();
      std::ofstream lineage(folder / "resume_source.txt");
      lineage << std::filesystem::absolute(source).string() << '\n';
      std::cout << "Resuming independently generated cold checkpoint: "
                << model(x) << std::endl;
    }
    bool final_converged = false;
    int evaluations = 0;
    // Small transformed gradients alone cannot certify saturated coordinates.
    // Reject a stopping point if an explicit inward probe lowers the objective.
    auto stationary = [&](const auto &checked) {
      if (!std::isfinite(checked.objective_value_m) ||
          !std::isfinite(checked.gradient_norm_m) ||
          checked.gradient_norm_m > 1e-4)
        return false;
      for (size_t i = 0; i < x.size(); ++i)
        if (std::abs(x[i]) > 16.) {
          auto inward = x;
          inward[i] = std::copysign(9.21024036697585, x[i]);
          if (model(inward) < checked.objective_value_m - 1e-8)
            return false;
        }
      return true;
    };
    for (int phase = resume.empty() ? (initial_valid ? 1 : 0) : 5; phase <= 4;
         ++phase) {
      swordfish::ControlStartPhase partial(model, x, std::max(1, phase),
                                           phase == 0);
      quadra::ParameterVector parameters;
      for (auto i : partial.active)
        parameters.add(quadra::Parameter(model.free[i].name, x[i]));
      std::vector<double> phase_start;
      for (auto i : partial.active)
        phase_start.push_back(x[i]);
      const double before = partial(phase_start);
      std::cout << "Activating phase " << phase << ", coordinates "
                << partial.active.size() << std::endl;
      const auto fit = quadra::optimize_lbfgs(
          partial, parameters, quadra::default_laplace_options(),
          phase < 3 ? 150 : (phase == 3 ? 150 : std::min(limit, 250)),
          phase < 4 ? 1e-3 : 1e-4);
      x = partial.expand(fit.par);
      evaluations += fit.iterations;
      for (const auto &r : fit.trajectory)
        trace << phase << ',' << r.evaluation << ',' << r.objective << ','
              << r.gradient_norm << '\n';
      trace.flush();
      stages << phase << ',' << partial.active.size() << ',' << before << ','
             << fit.value << ',' << fit.grad_norm << ',' << fit.converged << ','
             << fit.iterations << '\n';
      stages.flush();
      auto checkpoint = model.decode(x);
      save_parameters(
          folder / ("phase_" + std::to_string(phase) + "_parameters.csv"),
          checkpoint.first, checkpoint.second);
      final_converged = fit.converged;
    }
    // Learn scaling from this cold fit only. These eigenvalue adjustments
    // define an optimizer metric, not an inferential covariance or a changed
    // likelihood.
    for (int pass = 0; pass < 4; ++pass) {
      // A saturated logit can conceal a feasible inward descent direction.
      // Restore numerical sensitivity before learning a metric; this changes
      // only the start.
      std::ofstream recovery(
          folder / ("boundary_recovery_" + std::to_string(pass) + ".csv"));
      recovery << "parameter,old_logit,new_logit\n" << std::setprecision(17);
      for (size_t i = 0; i < x.size(); ++i)
        if (std::abs(x[i]) > 16.) {
          const double inward = std::copysign(9.21024036697585, x[i]);
          recovery << model.free[i].name << ',' << x[i] << ',' << inward
                   << '\n';
          x[i] = inward;
        }
      recovery.close();
      const int n = x.size();
      Eigen::MatrixXd hessian(n, n);
      for (int j = 0; j < n; ++j) {
        auto plus = x, minus = x;
        plus[j] += 1e-4;
        minus[j] -= 1e-4;
        const auto gp = quadra::evaluate_gradient(model, plus),
                   gm = quadra::evaluate_gradient(model, minus);
        for (int i = 0; i < n; ++i)
          hessian(i, j) = (gp.gradient_m[i] - gm.gradient_m[i]) / 2e-4;
        if (j % 14 == 0)
          std::cout << "Learning cold-fit scaling column " << j << '/' << n
                    << std::endl;
      }
      const Eigen::MatrixXd symmetric = (hessian + hessian.transpose()) * .5;
      Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> solver(symmetric);
      if (solver.info() != Eigen::Success || !solver.eigenvalues().allFinite())
        throw std::runtime_error("Cold curvature eigensolver failure");
      const double floor =
          std::max(1e-5, solver.eigenvalues().cwiseAbs().maxCoeff() * 1e-9);
      const Eigen::MatrixXd factor =
          solver.eigenvectors() * solver.eigenvalues()
                                      .cwiseAbs()
                                      .cwiseMax(floor)
                                      .cwiseSqrt()
                                      .cwiseInverse()
                                      .asDiagonal();
      std::ofstream metric(
          folder / ("optimizer_metric_" + std::to_string(pass) + ".json"));
      metric << std::setprecision(17)
             << "{\"source\":\"AD gradient differences at cold-fit "
                "checkpoint\",\"step\":0.0001,\"eigenvalue_floor\":"
             << floor << ",\"minimum_raw_eigenvalue\":"
             << solver.eigenvalues().minCoeff()
             << ",\"not_inferential_covariance\":true}\n";
      ColdMetricModel scaled(model, x, factor);
      quadra::ParameterVector scaled_parameters;
      for (int i = 0; i < n; ++i)
        scaled_parameters.add(quadra::Parameter(model.free[i].name, 0.));
      const auto refined = quadra::optimize_lbfgs(
          scaled, scaled_parameters, quadra::default_laplace_options(),
          std::min(limit, 200), 1e-7);
      x = scaled.expand(refined.par);
      evaluations += refined.iterations;
      final_converged = refined.converged;
      for (const auto &r : refined.trajectory)
        trace << (5 + pass) << ',' << r.evaluation << ',' << r.objective << ','
              << r.gradient_norm << '\n';
      stages << (5 + pass) << ',' << n << ',' << model(scaled.center) << ','
             << refined.value << ',' << refined.grad_norm << ','
             << refined.converged << ',' << refined.iterations << '\n';
      trace.flush();
      stages.flush();
      metric.flush();
      const auto current = model.decode(x);
      save_parameters(
          folder / ("refinement_" + std::to_string(pass) + "_parameters.csv"),
          current.first, current.second);
      const auto checked = quadra::evaluate_gradient(model, x);
      final_converged = stationary(checked);
      std::cout << "Refinement pass " << pass << ", original gradient "
                << checked.gradient_norm_m << std::endl;
      if (final_converged)
        break;
    }
    const auto gradient = quadra::evaluate_gradient(model, x);
    final_converged = stationary(gradient);
    auto fitted = model.decode(x);
    const auto likelihood = model.likelihood(x);
    for (const auto &row : likelihood.indices)
      fitted.first.at("LnQ_fleet_" + std::to_string(row.fleet)) = row.log_q;
    save_parameters(folder / "fitted_parameters.csv", fitted.first,
                    fitted.second);
    swordfish::write_fitted_outputs(data, fitted.first, fitted.second, folder);
    std::ofstream summary(folder / "summary.json");
    summary << std::setprecision(17) << "{\"mode\":\"" << mode
            << "\",\"initialization\":\"control native initial values; all "
               "historical deviations "
               "zero\",\"reads_SS3_estimates\":false,\"uses_SS3_covariance\":"
               "false,\"fixed_male_L1\":"
            << p.at("L_at_Amin_Mal_GP_1") << ",\"initial_objective\":"
            << (initial_valid ? std::to_string(initial) : "null")
            << ",\"literal_start_evaluable\":"
            << (initial_valid ? "true" : "false")
            << ",\"catch_feasibility_stage\":"
            << (initial_valid ? "false" : "true")
            << ",\"objective\":" << gradient.objective_value_m
            << ",\"gradient_norm\":" << gradient.gradient_norm_m
            << ",\"converged\":" << (final_converged ? "true" : "false")
            << ",\"free_parameters\":" << x.size()
            << ",\"optimizer_iterations\":" << evaluations << "}\n";
    if (!trace || !stages || !summary)
      throw std::runtime_error("Control-start output write failure");
    std::cout << "Final objective " << gradient.objective_value_m
              << ", full gradient " << gradient.gradient_norm_m
              << ", converged " << final_converged << std::endl;
    return final_converged ? 0 : 1;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 2;
  }
}
