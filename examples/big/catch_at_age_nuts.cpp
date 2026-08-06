#include "../../include/quadra/sampling.hpp"
#include "catch_at_age_shared.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

struct NoncenteredCatchAtAgePosterior {
  example::CatchAtAgeLaplaceModel *model;

  template <class T> T operator()(const std::vector<T> &q) const {
    const int years = model->data.n_years;
    std::vector<T> parameters(static_cast<std::size_t>(10 + years));
    for (int j = 0; j < 10; ++j) parameters[static_cast<std::size_t>(j)] = q[j];

    const T sigma_rec = example::positive_floor(q[8], 0.05);
    for (int y = 0; y < years; ++y) {
      parameters[static_cast<std::size_t>(10 + y)] =
          sigma_rec * q[static_cast<std::size_t>(10 + y)];
    }
    // model(parameters) is a density with respect to centered recruitment
    // deviations. Transforming u = sigma_rec * z adds this Jacobian.
    return -(*model)(parameters) + T(years) * log(sigma_rec);
  }
};

int main() {
  example::CatchAtAgeLaplaceModel model;
  // Keep the demonstration quick while exercising the same age-structured
  // process and robust Dirichlet-multinomial likelihood as the large fit.
  model.data.n_years = 6;

  std::vector<double> initial{
      std::log(900.0), std::log(0.25), std::log(0.15), std::log(0.18), 0.0,
      std::log(1.25),  std::log(0.20), std::log(0.18), std::log(0.35),
      std::log(40.0)};
  initial.resize(static_cast<std::size_t>(10 + model.data.n_years), 0.0);

  NoncenteredCatchAtAgePosterior posterior{&model};
  quadra::sampling::ReusableAdLogDensity<NoncenteredCatchAtAgePosterior>
      replay_check(posterior, initial);
  std::vector<double> replay_probe = initial;
  for (std::size_t i = 0; i < replay_probe.size(); ++i)
    replay_probe[i] += 1.0e-3 * static_cast<double>(i + 1);
  const auto replayed = replay_check.evaluate(replay_probe);
  const auto fresh =
      quadra::sampling::evaluate_ad_log_density(posterior, replay_probe);
  double max_gradient_difference = 0.0;
  for (std::size_t i = 0; i < replay_probe.size(); ++i)
    max_gradient_difference =
        std::max(max_gradient_difference,
                 std::abs(replayed.gradient[i] - fresh.gradient[i]));
  if (std::abs(replayed.log_density - fresh.log_density) > 1.0e-10 ||
      max_gradient_difference > 1.0e-10 ||
      replay_check.unsupported_replay_vertex_count() != 0)
    throw std::runtime_error("catch-at-age AD replay validation failed");
  quadra::sampling::NutsOptions options;
  options.warmup = 500;
  options.samples = 500;
  options.max_tree_depth = 10;
  options.target_acceptance = 0.80;
  options.seed = 20260806;
  std::vector<std::vector<double>> initial_states(4, initial);
  for (std::size_t chain = 0; chain < initial_states.size(); ++chain) {
    const double offset = 0.01 * (static_cast<double>(chain) - 1.5);
    for (std::size_t i = 0; i < initial_states[chain].size(); ++i)
      initial_states[chain][i] += offset / static_cast<double>(i + 1);
  }
  const auto fit = quadra::sampling::sample_nuts_chains(
      [&model](std::size_t) {
        return NoncenteredCatchAtAgePosterior{&model};
      },
      initial_states, options, true);

  double mean_log_m = 0.0;
  int total_draws = 0;
  int divergences = 0;
  int depth_hits = 0;
  for (const auto &chain : fit.chains) {
    for (const auto &draw : chain.draws) mean_log_m += draw[1];
    total_draws += static_cast<int>(chain.draws.size());
    divergences += chain.diagnostics.divergences;
    depth_hits += chain.diagnostics.max_depth_hits;
  }
  mean_log_m /= total_draws;
  const auto worst_rhat = std::max_element(fit.diagnostics.split_rhat.begin(),
                                           fit.diagnostics.split_rhat.end());
  const auto worst_bulk = std::min_element(fit.diagnostics.bulk_ess.begin(),
                                           fit.diagnostics.bulk_ess.end());
  const auto worst_tail = std::min_element(fit.diagnostics.tail_ess.begin(),
                                           fit.diagnostics.tail_ess.end());
  const std::size_t worst_rhat_parameter =
      static_cast<std::size_t>(worst_rhat - fit.diagnostics.split_rhat.begin());
  std::vector<std::string> parameter_names{
      "log_R0",          "log_M",        "log_F",       "log_q",
      "logit_sel50",     "log_sel_slope", "log_sigma_I", "log_sigma_C",
      "log_sigma_rec",   "log_comp_concentration"};
  for (int year = 0; year < model.data.n_years; ++year)
    parameter_names.push_back("z_rec[" + std::to_string(year) + "]");
  std::cout << "Quadra non-centered joint AD-NUTS catch-at-age example\n"
            << "chains: " << fit.chains.size() << "\n"
            << "draws: " << total_draws << "\n"
            << "mean log_M: " << mean_log_m << "\n"
            << "log_M split R-hat: " << fit.diagnostics.split_rhat[1] << "\n"
            << "log_M bulk ESS: " << fit.diagnostics.bulk_ess[1] << "\n"
            << "log_M tail ESS: " << fit.diagnostics.tail_ess[1] << "\n"
            << "maximum split R-hat: " << *worst_rhat << " ("
            << parameter_names[worst_rhat_parameter] << ")\n"
            << "minimum bulk ESS: " << *worst_bulk << "\n"
            << "minimum tail ESS: " << *worst_tail << "\n"
            << "divergences: " << divergences << "\n"
            << "max-depth hits: " << depth_hits << "\n";
  for (std::size_t chain = 0; chain < fit.chains.size(); ++chain) {
    const auto &diagnostics = fit.chains[chain].diagnostics;
    std::cout << "chain " << chain + 1
              << ": acceptance=" << diagnostics.mean_acceptance
              << " step_size=" << diagnostics.step_size
              << " divergences=" << diagnostics.divergences
              << " depth_hits=" << diagnostics.max_depth_hits
              << " BFMI=" << diagnostics.energy_bfmi << "\n";
  }
  return divergences > 20 || *worst_rhat > 1.2
             ? 1
             : 0;
}
