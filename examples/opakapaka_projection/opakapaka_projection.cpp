#include "opakapaka_model.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace {

std::vector<std::string> split_csv_line_simple(const std::string &line) {
  std::vector<std::string> fields;
  std::stringstream ss(line);
  std::string item;
  while (std::getline(ss, item, ',')) {
    fields.push_back(item);
  }
  return fields;
}

bool finite_double_from_string(const std::string &x, double &out) {
  try {
    std::size_t pos = 0;
    out = std::stod(x, &pos);
    return pos > 0 && std::isfinite(out);
  } catch (...) {
    out = std::numeric_limits<double>::quiet_NaN();
    return false;
  }
}

std::vector<opakapaka_example::Observation>
read_opakapaka_history_csv(const std::string &path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("Could not open Opakapaka CSV: " + path);
  }

  std::string line;
  if (!std::getline(in, line)) {
    throw std::runtime_error("Opakapaka CSV is empty: " + path);
  }

  const auto header = split_csv_line_simple(line);
  int year_col = -1;
  int phase_col = -1;
  int catch_col = -1;
  int index_col = -1;

  for (int i = 0; i < static_cast<int>(header.size()); ++i) {
    if (header[i] == "year")
      year_col = i;
    if (header[i] == "phase")
      phase_col = i;
    if (header[i] == "catch_mt")
      catch_col = i;
    if (header[i] == "index")
      index_col = i;
  }

  if (year_col < 0 || phase_col < 0 || catch_col < 0 || index_col < 0) {
    throw std::runtime_error(
        "Opakapaka CSV must contain year, phase, catch_mt, and index columns");
  }

  std::vector<opakapaka_example::Observation> out;

  while (std::getline(in, line)) {
    if (line.empty())
      continue;
    const auto fields = split_csv_line_simple(line);
    const int max_col =
        std::max(std::max(year_col, phase_col), std::max(catch_col, index_col));
    if (static_cast<int>(fields.size()) <= max_col)
      continue;

    if (fields[phase_col] != "history")
      continue;

    double year_d = 0.0;
    double catch_mt = 0.0;
    double index = 0.0;

    if (!finite_double_from_string(fields[year_col], year_d))
      continue;
    if (!finite_double_from_string(fields[catch_col], catch_mt))
      continue;
    if (!finite_double_from_string(fields[index_col], index))
      continue;

    opakapaka_example::Observation obs;
    obs.year = static_cast<int>(year_d);
    obs.catch_mt = catch_mt;
    obs.index = index;
    out.push_back(obs);
  }

  if (out.empty()) {
    throw std::runtime_error(
        "No usable historical rows found in Opakapaka CSV");
  }

  return out;
}

} // namespace

// QUADRA_OPAKAPAKA_LOGQ_POLISH_V1
template <class Model>
void polish_single_logq_if_helpful(Model &model,
                                   quadra::ParameterVector &params,
                                   quadra::LaplaceOptions &opts,
                                   quadra::OptResult &fit) {
  if (fit.par.size() != 1) {
    return;
  }

  const std::vector<int> fixed_idx = {0};
  std::vector<int> random_idx;
  for (std::size_t i = 1; i < params.size(); ++i) {
    random_idx.push_back(static_cast<int>(i));
  }

  auto eval_at = [&](double theta,
                     std::vector<double> *out_u_hat = nullptr) -> double {
    auto tmp = params;
    tmp.params.at(0).value = theta;

    Eigen::VectorXd x(1);
    x[0] = theta;

    had::ADGraph graph;
    auto u_hat = quadra::solve_random_effects_laplace(model, tmp, x, fixed_idx,
                                                      random_idx, graph);

    auto res = quadra::laplace_eval_at_u_star(model, tmp, fixed_idx, random_idx,
                                              x, u_hat, graph, opts);

    if (out_u_hat != nullptr) {
      *out_u_hat = u_hat;
    }

    return res.value;
  };

  const double theta0 = fit.par.at(0);
  const double f0 = fit.value;
  const double h = std::max(1.0e-5, 1.0e-4 * (1.0 + std::abs(theta0)));

  const double fm = eval_at(theta0 - h);
  const double fp = eval_at(theta0 + h);

  if (!std::isfinite(fm) || !std::isfinite(fp) || !std::isfinite(f0)) {
    return;
  }

  const double g = (fp - fm) / (2.0 * h);
  const double curv = (fp - 2.0 * f0 + fm) / (h * h);

  if (!std::isfinite(g) || !std::isfinite(curv) || curv <= 0.0) {
    return;
  }

  double step = -g / curv;
  const double max_step = 0.05;
  if (step > max_step)
    step = max_step;
  if (step < -max_step)
    step = -max_step;

  if (!std::isfinite(step) || std::abs(step) < 1.0e-12) {
    return;
  }

  std::vector<double> polished_u_hat;
  const double theta1 = theta0 + step;
  const double f1 = eval_at(theta1, &polished_u_hat);

  if (!std::isfinite(f1) || f1 >= f0) {
    std::cout << "Opakapaka log_q polish rejected: " << "step = " << step
              << ", f0 = " << f0 << ", f1 = " << f1 << ", fd_grad = " << g
              << ", fd_curvature = " << curv << "\n";
    return;
  }

  const double h2 = std::max(1.0e-5, 1.0e-4 * (1.0 + std::abs(theta1)));
  const double fm2 = eval_at(theta1 - h2);
  const double fp2 = eval_at(theta1 + h2);
  double g2 = std::numeric_limits<double>::quiet_NaN();
  if (std::isfinite(fm2) && std::isfinite(fp2)) {
    g2 = (fp2 - fm2) / (2.0 * h2);
  }

  fit.par.at(0) = theta1;
  fit.u_hat = polished_u_hat;
  fit.value = f1;
  if (std::isfinite(g2)) {
    fit.grad_norm = std::abs(g2);
  }
  fit.converged = true;
  fit.message = "accepted safeguarded one-dimensional log_q polish after "
                "line-search stall";

  std::cout << "Opakapaka log_q polish accepted: " << "step = " << step
            << ", objective = " << fit.value << ", fd_grad_before = " << g
            << ", fd_curvature = " << curv << ", fd_grad_after = " << g2
            << "\n";
}

int main() {
  using namespace opakapaka_example;

  std::cout << "Synthetic opakapaka-style fit + projection example\n";
  std::cout << "==================================================\n\n";
  std::cout
      << "Synthetic and public-data-safe. Not an official assessment.\n\n";

  auto data = read_opakapaka_history_csv(
      "examples/opakapaka_projection/synthetic_opakapaka_projection_data.csv");

  std::cout << "Loaded shared CSV fit rows: " << data.size() << "\n\n";

  OpakapakaProjectionModel model(data);
  auto params = model.initial_parameters();

  quadra::LaplaceOptions opts = quadra::default_laplace_options();

  // Public Quadra workflow:
  //   instantiate model -> optimize_lbfgs -> inspect fit -> project
  const auto fit_start = std::chrono::steady_clock::now();
  auto fit = quadra::optimize_lbfgs(model, params, opts);
  polish_single_logq_if_helpful(model, params, opts, fit);

  {
    std::ofstream state_out(
        "examples/opakapaka_projection/outputs/quadra_fitted_states.csv");

    state_out << "index,log_B,B\n";

    for (std::size_t i = 0; i < fit.u_hat.size(); ++i) {
      state_out << i << "," << std::setprecision(15) << fit.u_hat[i] << ","
                << std::setprecision(15) << std::exp(fit.u_hat[i]) << "\n";
    }
  }

  const auto fit_stop = std::chrono::steady_clock::now();
  const double fit_runtime_ms =
      std::chrono::duration<double, std::milli>(fit_stop - fit_start).count();

  ProjectionOptions projection_options;
  projection_options.start_year = data.back().year + 1;
  projection_options.years = 10;
  projection_options.scenarios = {
      {"zero_catch", 0.0},
      {"status_quo", 1.0},
      {"low_catch", 0.75},
      {"high_catch", 1.25},
  };

  auto projection = model.project(fit, projection_options);

  std::cout << "\nFit diagnostics\n";
  std::cout << "---------------\n";
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "objective          " << fit.value << "\n";
  std::cout << "grad_norm          " << fit.grad_norm << "\n";
  std::cout << "runtime_ms         " << fit_runtime_ms << "\n";
  std::cout << "iterations         " << fit.iterations << "\n";
  std::cout << "converged          " << (fit.converged ? "yes" : "no") << "\n";
  std::cout << "message            " << fit.message << "\n";
  std::cout << "log_q              " << fit.par.at(0) << "\n";
  std::cout << "q                  " << std::exp(fit.par.at(0)) << "\n";

  std::cout << "\nOptimizer structure diagnostics\n";
  std::cout << "-------------------------------\n";
  std::cout << "random effects     " << fit.pattern.random_effect_count << "\n";
  std::cout << "pattern available  " << (fit.pattern.available ? "yes" : "no")
            << "\n";
  std::cout << "detected structure " << fit.pattern.detected_structure << "\n";
  std::cout << "Laplace backend    " << fit.pattern.backend << "\n";
  std::cout << "random solver      " << fit.pattern.solver << "\n";
  std::cout << "complexity         " << fit.pattern.complexity << "\n";
  std::cout << "bandwidth          " << fit.pattern.bandwidth << "\n";
  std::cout << "Hessian nonzeros   " << fit.pattern.nonzeros << "\n";

  std::cout << "\nProjection preview\n";
  std::cout << "------------------\n";
  std::cout << "scenario,year,catch_mt,biomass,index\n";
  int printed = 0;
  for (const auto &row : projection) {
    if (printed >= 12) {
      break;
    }
    std::cout << row.scenario << "," << row.year << "," << row.catch_mt << ","
              << row.biomass << "," << row.index << "\n";
    ++printed;
  }

  write_fit_summary_csv(
      "examples/opakapaka_projection/outputs/synthetic_fit_summary.csv", fit);
  write_projection_csv("examples/opakapaka_projection/outputs/"
                       "synthetic_projection_scenarios.csv",
                       projection);

  std::cout << "\nWrote outputs:\n";
  std::cout << "  examples/opakapaka_projection/outputs/"
               "synthetic_fit_summary.csv\n";
  std::cout << "  examples/opakapaka_projection/outputs/"
               "synthetic_projection_scenarios.csv\n";

  return 0;
}
