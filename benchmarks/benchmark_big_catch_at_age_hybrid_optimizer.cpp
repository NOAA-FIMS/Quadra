#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <sys/resource.h>
#include <vector>

#include "../examples/big/catch_at_age_shared.hpp"
#include "../include/quadra/stats.hpp"

namespace {
double peak_rss_mb() {
  rusage usage{};
  getrusage(RUSAGE_SELF, &usage);
#ifdef __APPLE__
  return static_cast<double>(usage.ru_maxrss) / (1024.0 * 1024.0);
#else
  return static_cast<double>(usage.ru_maxrss) / 1024.0;
#endif
}
} // namespace

int main(int argc, char **argv) {
  using Clock = std::chrono::steady_clock;
  const std::string mode = argc > 1 ? argv[1] : "hybrid";
  const double switch_threshold = argc > 2 ? std::strtod(argv[2], nullptr) : 0.1;
  if (mode != "hybrid" && mode != "exact") {
    std::cerr << "usage: benchmark_big_catch_at_age_hybrid_optimizer "
                 "[hybrid|exact] [switch_threshold]\n";
    return 2;
  }

  example::CatchAtAgeLaplaceModel model;
  quadra::ParameterSet parameters;
  const std::vector<double> initial = {
      std::log(900.0), std::log(0.25), std::log(0.15), std::log(0.18), 0.0,
      std::log(1.25), std::log(0.20), std::log(0.15), std::log(0.35),
      std::log(40.0)};
  const std::vector<const char *> names = {
      "log_R0",          "log_M",           "log_q",
      "log_Fbar",        "sel50_raw",       "log_sel_slope",
      "log_sigma_index", "log_sigma_catch", "log_sigma_rec",
      "log_comp_concentration"};
  for (std::size_t i = 0; i < initial.size(); ++i)
    parameters.add(names[i], initial[i], quadra::ParameterTransform::Identity,
                   false);
  std::vector<double> random(static_cast<std::size_t>(model.data.n_years), 0.0);
  for (int year = 0; year < model.data.n_years; ++year)
    parameters.add("rec_dev_" + std::to_string(year + 1), 0.0,
                   quadra::ParameterTransform::Identity, true);

  quadra::LaplaceObjectiveOptions objective_options;
  objective_options.include_constant_m = true;
  objective_options.newton_m.gradient_tolerance_m = 1.0e-6;
  quadra::laplace::ExactLaplaceGradientEngineOptions exact_engine;
  exact_engine.discover_active_directions = false;
  exact_engine.stream_dense_hdot_trace = true;
  exact_engine.dense_hdot_bandwidth = -1;
  quadra::stats::ExactLaplaceEvaluator<example::CatchAtAgeLaplaceModel> exact(
      model, initial, random, parameters, objective_options, exact_engine);

  quadra::stats::LaplaceOptimizerOptions optimizer_options;
  optimizer_options.max_iterations = 100;
  optimizer_options.memory = 20;
  optimizer_options.gradient_tolerance = 1.0e-4;
  optimizer_options.maximum_direction_norm = 1.0;

  const auto start = Clock::now();
  quadra::stats::LaplaceOptimizerResult result;
  if (mode == "exact") {
    result = quadra::stats::optimize_laplace(exact, initial, optimizer_options);
  } else {
    result = quadra::stats::optimize_laplace_hybrid(
        exact, initial, 5, switch_threshold, optimizer_options);
  }
  const double wall_ms =
      std::chrono::duration<double, std::milli>(Clock::now() - start).count();

  std::cout << "mode,switch_threshold,converged,iterations,approximate_evaluations,"
               "exact_evaluations,switch_iteration,objective,gradient_norm,"
               "wall_ms,peak_rss_mb,message\n";
  std::cout << mode << ',' << switch_threshold << ','
            << (result.converged ? 1 : 0) << ','
            << result.iterations << ',' << result.approximate_evaluations << ','
            << result.exact_evaluations << ',' << result.exact_switch_iteration
            << ',' << std::setprecision(15) << result.objective << ','
            << result.gradient_norm << ',' << std::fixed << std::setprecision(3)
            << wall_ms << ',' << peak_rss_mb() << ',' << result.message << '\n';
}
