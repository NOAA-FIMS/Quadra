#include "../include/quadra/sampling.hpp"

#include <cmath>
#include <iostream>
#include <sstream>
#include <vector>

DECLARE_ADGRAPH();

struct CorrelatedGaussian {
  template <class T> T operator()(const std::vector<T> &q) const {
    const T x = q[0] - T(1.0);
    const T y = q[1] + T(0.5);
    return -T(0.5) * (T(1.3333333333333333) * x * x -
                      T(1.3333333333333333) * x * y +
                      T(1.3333333333333333) * y * y);
  }
};

struct AnisotropicGaussian {
  template <class T> T operator()(const std::vector<T> &q) const {
    return -T(0.5) * (q[0] * q[0] / T(9.0) + q[1] * q[1] / T(0.04));
  }
};

struct HighlyCorrelatedGaussian {
  template <class T> T operator()(const std::vector<T> &q) const {
    constexpr double rho = 0.99;
    constexpr double scale = 1.0 / (1.0 - rho * rho);
    return -T(0.5 * scale) *
           (q[0] * q[0] - T(2.0 * rho) * q[0] * q[1] + q[1] * q[1]);
  }
};

struct StandardGaussian {
  template <class T> T operator()(const std::vector<T> &q) const {
    T log_density = T(0.0);
    for (const T &value : q) log_density -= T(0.5) * value * value;
    return log_density;
  }
};

struct BoundedTarget {
  template <class T> T operator()(const std::vector<T> &q) const {
    if (quadra::value_of(q[0]) > 2.0)
      throw std::domain_error("outside support");
    return -T(0.5) * q[0] * q[0];
  }
};

struct SquareReplayTarget {
  template <class T> T operator()(const std::vector<T> &q) const {
    using had::square;
    return -T(0.5) * square(q[0]);
  }
};

int main() {
  CorrelatedGaussian target;
  quadra::sampling::NutsOptions options;
  options.warmup = 500;
  options.samples = 1000;
  options.max_tree_depth = 8;
  options.seed = 20260806;
  const auto result =
      quadra::sampling::sample_nuts(target, std::vector<double>{0.0, 0.0},
                                   options);
  AnisotropicGaussian anisotropic_target;
  quadra::sampling::NutsOptions anisotropic_options = options;
  anisotropic_options.seed = 20260807;
  const auto anisotropic = quadra::sampling::sample_nuts(
      anisotropic_target, std::vector<double>{0.5, 0.5},
      anisotropic_options);
  HighlyCorrelatedGaussian highly_correlated_target;
  quadra::sampling::NutsOptions correlated_diagonal_options = options;
  correlated_diagonal_options.seed = 20260808;
  const auto correlated_diagonal = quadra::sampling::sample_nuts(
      highly_correlated_target, std::vector<double>{0.0, 0.0},
      correlated_diagonal_options);
  quadra::sampling::NutsOptions correlated_dense_options =
      correlated_diagonal_options;
  correlated_dense_options.adapt_dense_mass = true;
  const auto correlated_dense = quadra::sampling::sample_nuts(
      highly_correlated_target, std::vector<double>{0.0, 0.0},
      correlated_dense_options);
  StandardGaussian standard_target;
  quadra::sampling::NutsOptions short_dense_options;
  short_dense_options.warmup = 20;
  short_dense_options.samples = 25;
  short_dense_options.adapt_dense_mass = true;
  short_dense_options.seed = 20260809;
  const auto short_dense = quadra::sampling::sample_nuts(
      standard_target, std::vector<double>(12, 0.1), short_dense_options);
  quadra::sampling::NutsOptions minimal_warmup_options;
  minimal_warmup_options.warmup = 1;
  minimal_warmup_options.samples = 5;
  minimal_warmup_options.adapt_dense_mass = true;
  const auto minimal_warmup = quadra::sampling::sample_nuts(
      standard_target, std::vector<double>(2, 0.1), minimal_warmup_options);
  quadra::sampling::detail::EuclideanMetric rejected_metric(2);
  const bool accepted_indefinite = rejected_metric.set_inverse_mass(
      {1.0, 2.0, 2.0, 1.0}, true);
  quadra::sampling::NutsOptions multi_options = options;
  multi_options.warmup = 300;
  multi_options.samples = 500;
  const std::vector<std::vector<double>> initial_states{
      {-1.0, -1.0}, {0.0, 0.0}, {2.0, -1.0}, {1.0, 1.0}};
  const auto multi = quadra::sampling::sample_nuts_chains(
      [](std::size_t) { return CorrelatedGaussian{}; }, initial_states,
      multi_options, true);
  quadra::sampling::NutsHealthThresholds health_thresholds;
  health_thresholds.max_rhat = 1.05;
  health_thresholds.max_divergences = 5;
  health_thresholds.max_depth_hits = 5;
  const auto health =
      quadra::sampling::assess_nuts_health(multi, health_thresholds);
  quadra::sampling::NutsWorkflowOptions workflow_options;
  workflow_options.sampler = multi_options;
  workflow_options.sampler.warmup = 150;
  workflow_options.sampler.samples = 150;
  workflow_options.health.max_rhat = 1.1;
  workflow_options.health.min_bulk_ess = 20.0;
  workflow_options.health.min_tail_ess = 20.0;
  workflow_options.health.max_divergences = 10;
  workflow_options.health.max_depth_hits = 10;
  workflow_options.initialization_seed = 20260811;
  const auto workflow = quadra::sampling::run_nuts_workflow(
      [](std::size_t) { return CorrelatedGaussian{}; }, {1.0, -0.5},
      {"x", "y"}, workflow_options);
  quadra::sampling::PosteriorSimulationOptions simulation_options;
  simulation_options.thin = 3;
  simulation_options.max_draws = 7;
  simulation_options.seed = 20260812;
  auto simulator = [](const std::vector<double> &draw, std::mt19937_64 &rng) {
    std::normal_distribution<double> normal;
    return std::vector<double>{draw[0] + normal(rng), draw[1]};
  };
  const auto simulation = quadra::sampling::simulate_posterior(
      workflow, {"replicated_x,raw", "latent_y"}, simulator,
      simulation_options);
  const auto repeated_simulation = quadra::sampling::simulate_posterior(
      workflow, {"replicated_x,raw", "latent_y"}, simulator,
      simulation_options);
  std::ostringstream posterior_csv;
  std::ostringstream parameter_csv;
  std::ostringstream chain_csv;
  std::ostringstream predictive_csv;
  std::ostringstream summary_csv;
  quadra::sampling::write_posterior_draws_csv(posterior_csv, workflow);
  quadra::sampling::write_parameter_diagnostics_csv(parameter_csv, workflow);
  quadra::sampling::write_chain_diagnostics_csv(chain_csv, workflow);
  quadra::sampling::write_posterior_predictive_csv(predictive_csv, simulation);
  quadra::sampling::write_nuts_summary_csv(summary_csv, workflow);
  bool rejected_duplicate_names = false;
  try {
    quadra::sampling::validate_parameter_names({"x", "x"}, 2);
  } catch (const std::invalid_argument &) {
    rejected_duplicate_names = true;
  }
  std::vector<quadra::sampling::NutsResult> scale_mismatch(4);
  for (std::size_t chain = 0; chain < scale_mismatch.size(); ++chain) {
    const double scale = chain < 2 ? 1.0 : 3.0;
    for (int draw = 0; draw < 200; ++draw)
      scale_mismatch[chain].draws.push_back(
          {scale * std::sin(0.37 * draw + 0.11 * chain)});
  }
  const auto scale_diagnostics =
      quadra::sampling::compute_multi_chain_diagnostics(scale_mismatch);
  BoundedTarget bounded;
  quadra::sampling::NutsOptions bounded_options;
  bounded_options.warmup = 25;
  bounded_options.samples = 25;
  bounded_options.seed = 17;
  bounded_options.reuse_ad_graph = false;
  const auto bounded_result = quadra::sampling::sample_nuts(
      bounded, std::vector<double>{0.0}, bounded_options);
  SquareReplayTarget square_target;
  quadra::sampling::ReusableAdLogDensity<SquareReplayTarget> square_replay(
      square_target, {0.0});
  const auto square_replayed = square_replay.evaluate({0.3});
  const auto square_fresh =
      quadra::sampling::evaluate_ad_log_density(square_target, {0.3});
  std::vector<double> mean(2, 0.0);
  for (const auto &draw : result.draws) {
    mean[0] += draw[0];
    mean[1] += draw[1];
  }
  mean[0] /= result.draws.size();
  mean[1] /= result.draws.size();
  if (result.draws.size() != 1000 || std::abs(mean[0] - 1.0) > 0.12 ||
      std::abs(mean[1] + 0.5) > 0.12 ||
      result.diagnostics.divergences > 5 ||
      !(result.diagnostics.mean_acceptance > 0.6) ||
      result.diagnostics.mass_matrix_updates < 2 ||
      anisotropic.diagnostics.inverse_mass.size() != 2 ||
      anisotropic.diagnostics.inverse_mass[0] <
          20.0 * anisotropic.diagnostics.inverse_mass[1] ||
      anisotropic.diagnostics.max_depth_hits > 5 ||
      correlated_dense.diagnostics.dense_inverse_mass.size() != 2 ||
      correlated_dense.diagnostics.dense_inverse_mass[0][1] < 0.5 ||
      correlated_dense.diagnostics.divergences > 5 ||
      correlated_dense.diagnostics.leapfrog_steps >=
          correlated_diagonal.diagnostics.leapfrog_steps ||
      short_dense.draws.size() != 25 ||
      short_dense.diagnostics.mass_matrix_updates < 1 ||
      short_dense.diagnostics.dense_inverse_mass.size() != 12 ||
      minimal_warmup.draws.size() != 5 ||
      minimal_warmup.diagnostics.mass_matrix_updates != 0 ||
      accepted_indefinite || rejected_metric.dense ||
      rejected_metric.inverse_mass[0] != 1.0 ||
      rejected_metric.inverse_mass[3] != 1.0 ||
      !(result.diagnostics.energy_bfmi > 0.0) ||
      !health.passed ||
      !workflow.health.passed || workflow.parameter_count() != 2 ||
      workflow.total_draws() != 600 ||
      workflow.initial_states.size() != 4 || !rejected_duplicate_names ||
      simulation.draws.size() != 7 ||
      simulation.quantity_names.size() != 2 ||
      simulation.draws[1].iteration != 3 ||
      simulation.draws[0].values != repeated_simulation.draws[0].values ||
      posterior_csv.str().find("chain,iteration,parameter") != 0 ||
      parameter_csv.str().find("parameter,rhat,bulk_ess,tail_ess") != 0 ||
      chain_csv.str().find("mass_update_failures") == std::string::npos ||
      predictive_csv.str().find("\"replicated_x,raw\"") ==
          std::string::npos ||
      summary_csv.str().find("health,PASS") == std::string::npos ||
      multi.chains.size() != 4 || multi.diagnostics.split_rhat[0] > 1.05 ||
      multi.diagnostics.split_rhat[1] > 1.05 ||
      multi.diagnostics.bulk_ess[0] < 100.0 ||
      multi.diagnostics.bulk_ess[1] < 100.0 ||
      multi.diagnostics.tail_ess[0] < 50.0 ||
      multi.diagnostics.tail_ess[1] < 50.0 ||
      scale_diagnostics.split_rhat[0] < 1.1 ||
      bounded_result.draws.size() != 25 ||
      std::abs(square_replayed.log_density - square_fresh.log_density) >
          1.0e-14 ||
      std::abs(square_replayed.gradient[0] - square_fresh.gradient[0]) >
          1.0e-14 ||
      square_replay.unsupported_replay_vertex_count() != 0) {
    std::cerr << "FAIL: AD NUTS Gaussian recovery failed\n"
              << "mean=" << mean[0] << "," << mean[1]
              << " acceptance=" << result.diagnostics.mean_acceptance
              << " divergences=" << result.diagnostics.divergences << "\n";
    return 1;
  }
  std::cout << "PASS: AD NUTS Gaussian recovery\n"
            << "  mean: " << mean[0] << ", " << mean[1] << "\n"
            << "  acceptance: " << result.diagnostics.mean_acceptance << "\n"
            << "  step size: " << result.diagnostics.step_size << "\n";
  std::cout << "  split R-hat: " << multi.diagnostics.split_rhat[0] << ", "
            << multi.diagnostics.split_rhat[1] << "\n"
            << "  bulk ESS: " << multi.diagnostics.bulk_ess[0] << ", "
            << multi.diagnostics.bulk_ess[1] << "\n";
  return 0;
}
