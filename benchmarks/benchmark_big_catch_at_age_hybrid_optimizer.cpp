#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <sys/resource.h>
#include <vector>

#include "../core/inference/fixed_effect_covariance.hpp"
#include "../core/laplace/random_effect_hessian.hpp"
#include "../core/laplace/laplace_exact_directional_curvature.hpp"
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
  const int years = argc > 3 ? std::atoi(argv[3]) : 30;
  if (mode != "hybrid" && mode != "exact") {
    std::cerr << "usage: benchmark_big_catch_at_age_hybrid_optimizer "
                 "[hybrid|exact] [switch_threshold] [years]\n";
    return 2;
  }

  example::CatchAtAgeLaplaceModel model(years);
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

  std::cout << "mode,switch_threshold,years,converged,iterations,approximate_evaluations,"
               "exact_evaluations,switch_iteration,objective,gradient_norm,"
               "wall_ms,peak_rss_mb,message\n";
  std::cout << mode << ',' << switch_threshold << ',' << years << ','
            << (result.converged ? 1 : 0) << ','
            << result.iterations << ',' << result.approximate_evaluations << ','
            << result.exact_evaluations << ',' << result.exact_switch_iteration
            << ',' << std::setprecision(15) << result.objective << ','
            << result.gradient_norm << ',' << std::fixed << std::setprecision(3)
            << wall_ms << ',' << peak_rss_mb() << ',' << result.message << '\n';

  if (!result.converged)
    return 1;

  auto exact_gradient = [&](const std::vector<double> &fixed) {
    return exact.evaluate(fixed, result.random_mode).gradient;
  };
  const auto fd_start = Clock::now();
  const auto fd = quadra::estimate_fixed_effect_covariance_from_gradient(
      exact_gradient, result.fixed, 1.0e-4);
  const double fd_ms =
      std::chrono::duration<double, std::milli>(Clock::now() - fd_start)
          .count();

  const auto schur_start = Clock::now();
  quadra::RandomEffectHessianWorkspace<example::CatchAtAgeLaplaceModel>
      workspace(model, result.fixed, result.random_mode,
                quadra::partition_parameters(parameters), false);
  const auto blocks =
      workspace.EvaluateJointHessianBlocks(result.fixed, result.random_mode);
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> factorization;
  factorization.compute(blocks.random_hessian_m);
  const Eigen::MatrixXd profiled_joint_hessian =
      blocks.fixed_hessian_m - blocks.mixed_hessian_m.transpose() *
                                   factorization.solve(blocks.mixed_hessian_m);
  const auto schur = quadra::estimate_fixed_effect_covariance_from_hessian(
      profiled_joint_hessian);
  const double schur_ms =
      std::chrono::duration<double, std::milli>(Clock::now() - schur_start)
          .count();

  double hessian_relative_error = std::nan("");
  double se_max_relative_error = std::nan("");
  if (fd.success_m && schur.success_m) {
    hessian_relative_error =
        (fd.hessian_m - schur.hessian_m).norm() / fd.hessian_m.norm();
    se_max_relative_error = 0.0;
    for (Eigen::Index i = 0; i < fd.covariance_m.rows(); ++i) {
      const double fd_se = std::sqrt(fd.covariance_m(i, i));
      const double schur_se = std::sqrt(schur.covariance_m(i, i));
      se_max_relative_error =
          std::max(se_max_relative_error, std::abs(schur_se / fd_se - 1.0));
    }
  }
  std::cout << "covariance_method,years,success,wall_ms,gradient_evaluations,"
               "hessian_relative_difference,max_se_relative_difference\n";
  std::cout << std::setprecision(10);
  std::cout << "fd_exact_laplace," << years << ',' << (fd.success_m ? 1 : 0)
            << ',' << fd_ms << ',' << 2 * result.fixed.size() << ",0,0\n";
  std::cout << "ad_profiled_joint_schur," << years << ','
            << (schur.success_m ? 1 : 0) << ',' << schur_ms << ",0,"
            << hessian_relative_error << ',' << se_max_relative_error << '\n';

  const auto exact_hessian = quadra::laplace::exact_laplace_hessian_fourth_order(
      model, result.fixed, result.random_mode,
      quadra::partition_parameters(parameters), 4);
  const auto exact_covariance =
      quadra::estimate_fixed_effect_covariance_from_hessian(
          exact_hessian.hessian_m);
  const double exact_hessian_relative_error =
      (exact_hessian.hessian_m - fd.hessian_m).norm() / fd.hessian_m.norm();
  double exact_se_relative_error = 0.0;
  for (Eigen::Index i = 0; i < fd.covariance_m.rows(); ++i) {
    const double fd_se = std::sqrt(fd.covariance_m(i, i));
    const double exact_se = std::sqrt(exact_covariance.covariance_m(i, i));
    exact_se_relative_error = std::max(
        exact_se_relative_error, std::abs(exact_se / fd_se - 1.0));
  }
  std::cout << "exact_fourth_ad," << years << ','
            << (exact_covariance.success_m ? 1 : 0) << ','
            << exact_hessian.total_ms_m << ",0,"
            << exact_hessian_relative_error << ',' << exact_se_relative_error
            << '\n';
}
