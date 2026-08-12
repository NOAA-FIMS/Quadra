#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <sys/resource.h>
#include <vector>

#include "../examples/big/catch_at_age_shared.hpp"
#include "../include/quadra/stats.hpp"

namespace {

double peak_rss_mb() {
  rusage usage{};
  if (getrusage(RUSAGE_SELF, &usage) != 0)
    return std::nan("");
#ifdef __APPLE__
  return static_cast<double>(usage.ru_maxrss) / (1024.0 * 1024.0);
#else
  return static_cast<double>(usage.ru_maxrss) / 1024.0;
#endif
}

} // namespace

int main(int argc, char **argv) {
  using Clock = std::chrono::steady_clock;
  const int workers = argc > 1 ? std::atoi(argv[1]) : 1;
  const int repetitions = argc > 2 ? std::atoi(argv[2]) : 3;
  const bool stream_dense_trace = argc > 3 ? std::atoi(argv[3]) != 0 : true;
  if (workers <= 0 || repetitions <= 0) {
    std::cerr << "usage: benchmark_big_catch_at_age_hdot_workers "
                 "[workers] [repetitions]\n";
    return 2;
  }

  example::CatchAtAgeLaplaceModel model;
  quadra::ParameterSet parameters;
  const std::vector<double> fixed = {8.199973,  -0.369740, 0.005970,  -1.194726,
                                     0.293207,  0.350404,  -1.709454, -2.741071,
                                     -3.407494, 7.595303};
  const std::vector<const char *> names = {
      "log_R0",          "log_M",
      "log_q",           "log_Fbar",
      "sel50_raw",       "log_sel_slope",
      "log_sigma_index", "log_sigma_catch",
      "log_sigma_rec",   "log_comp_concentration"};
  for (std::size_t i = 0; i < fixed.size(); ++i) {
    parameters.add(names[i], fixed[i], quadra::ParameterTransform::Identity,
                   false);
  }
  std::vector<double> random(static_cast<std::size_t>(model.data.n_years), 0.0);
  for (int year = 0; year < model.data.n_years; ++year) {
    parameters.add("rec_dev_" + std::to_string(year + 1), 0.0,
                   quadra::ParameterTransform::Identity, true);
  }

  quadra::LaplaceObjectiveOptions objective_options;
  objective_options.include_constant_m = true;
  objective_options.newton_m.gradient_tolerance_m = 1.0e-6;
  quadra::laplace::ExactLaplaceGradientEngineOptions engine_options;
  engine_options.discover_active_directions = false;
  engine_options.hdot_workers = workers;
  engine_options.stream_dense_hdot_trace = stream_dense_trace;

  const auto construct_start = Clock::now();
  quadra::stats::ExactLaplaceEvaluator<example::CatchAtAgeLaplaceModel>
      evaluator(model, fixed, random, parameters, objective_options,
                engine_options);
  const double construct_ms =
      std::chrono::duration<double, std::milli>(Clock::now() - construct_start)
          .count();

  double total_ms = 0.0;
  double objective_ms = 0.0;
  double factorization_ms = 0.0;
  double sensitivity_ms = 0.0;
  double hdot_ms = 0.0;
  double trace_ms = 0.0;
  double gradient_checksum = 0.0;
  for (int repetition = 0; repetition < repetitions; ++repetition) {
    const auto result = evaluator.evaluate(fixed);
    if (!result.success) {
      std::cerr << "exact evaluation failed: " << result.objective.message_m
                << '\n';
      return 1;
    }
    total_ms += result.timings.total_ms;
    objective_ms += result.timings.objective_ms;
    factorization_ms += result.timings.factorization_ms;
    sensitivity_ms += result.timings.mode_sensitivity_ms;
    hdot_ms += result.timings.hdot_ms;
    trace_ms += result.timings.trace_ms;
    gradient_checksum = 0.0;
    for (double value : result.gradient)
      gradient_checksum += value;
  }

  std::cout
      << "stream_dense_trace,requested_workers,actual_workers,construct_ms,"
         "mean_objective_ms,mean_factorization_ms,mean_sensitivity_ms,"
         "mean_hdot_ms,mean_trace_ms,mean_total_ms,peak_rss_mb,"
         "gradient_checksum\n";
  std::cout << (stream_dense_trace ? 1 : 0) << ',' << workers << ','
            << evaluator.hdot_worker_count() << ',' << std::fixed
            << std::setprecision(3) << construct_ms << ','
            << objective_ms / repetitions << ','
            << factorization_ms / repetitions << ','
            << sensitivity_ms / repetitions << ',' << hdot_ms / repetitions
            << ',' << trace_ms / repetitions << ',' << total_ms / repetitions
            << ',' << peak_rss_mb() << ',' << std::setprecision(12)
            << gradient_checksum << '\n';
}
