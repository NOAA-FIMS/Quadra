#include "../core/optimizer.hpp"
#include "../include/quadra/sampling.hpp"

#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <sstream>
#include <string>
#include <vector>

DECLARE_ADGRAPH();

struct HierarchicalNormalModel {
  std::vector<std::vector<double>> observations;

  template <class T> T operator()(const std::vector<T> &q) const {
    const T mu = q[0];
    const T log_sigma = q[1];
    const T tau = T(0.4);
    const T sigma = exp(log_sigma);
    T nll = T(0.5) * mu * mu / T(4.0);
    nll += T(0.5) * (log_sigma - T(std::log(0.3))) *
           (log_sigma - T(std::log(0.3))) / T(0.49);
    for (std::size_t group = 0; group < observations.size(); ++group) {
      const T z = q[2 + group];
      const T group_mean = mu + tau * z;
      nll += T(0.5) * z * z;
      for (double observation : observations[group]) {
        const T residual = (T(observation) - group_mean) / sigma;
        nll += T(0.5) * residual * residual + log_sigma;
      }
    }
    return nll;
  }
};

struct HierarchicalPosterior {
  HierarchicalNormalModel *model;
  template <class T> T operator()(const std::vector<T> &q) const {
    return -(*model)(q);
  }
};

int main() {
  constexpr double true_mu = 1.25;
  constexpr double true_tau = 0.4;
  constexpr double true_sigma = 0.3;
  const std::vector<double> true_z{-1.0, 0.0, 0.8};
  std::mt19937_64 data_rng(20260814);
  std::normal_distribution<double> normal;
  HierarchicalNormalModel model;
  model.observations.assign(3, {});
  for (std::size_t group = 0; group < model.observations.size(); ++group)
    for (int replicate = 0; replicate < 8; ++replicate)
      model.observations[group].push_back(
          true_mu + true_tau * true_z[group] + true_sigma * normal(data_rng));

  quadra::ParameterVector parameters;
  parameters.add({"mu", 0.0, quadra::ParameterTransform::Identity, false});
  parameters.add({"log_sigma", std::log(0.3),
                  quadra::ParameterTransform::Identity, false});
  for (std::size_t group = 0; group < model.observations.size(); ++group)
    parameters.add({"z[" + std::to_string(group) + "]", 0.0,
                    quadra::ParameterTransform::Identity, false});
  const auto mode =
      quadra::optimize_lbfgs(model, parameters,
                             quadra::default_laplace_options(), 100, 1.0e-6);
  if (!mode.converged || mode.par.size() != 5)
    throw std::runtime_error("hierarchical mode fit failed");

  quadra::sampling::NutsWorkflowOptions workflow_options;
  workflow_options.sampler.warmup = 400;
  workflow_options.sampler.samples = 400;
  workflow_options.sampler.target_acceptance = 0.85;
  workflow_options.sampler.adapt_dense_mass = true;
  workflow_options.sampler.seed = 20260815;
  workflow_options.initialization_seed = 20260816;
  workflow_options.health.max_rhat = 1.05;
  workflow_options.health.min_bulk_ess = 100.0;
  workflow_options.health.min_tail_ess = 100.0;
  const auto posterior = quadra::sampling::run_nuts_workflow(
      [&model](std::size_t) { return HierarchicalPosterior{&model}; }, mode.par,
      {"mu", "log_sigma", "z[0]", "z[1]", "z[2]"},
      workflow_options);

  quadra::sampling::PosteriorSimulationOptions simulation_options;
  simulation_options.thin = 8;
  simulation_options.max_draws = 100;
  simulation_options.seed = 20260817;
  if (!posterior.health.passed)
    std::cerr << "sampling health: rhat=" << posterior.health.max_rhat
              << " bulk=" << posterior.health.min_bulk_ess
              << " tail=" << posterior.health.min_tail_ess
              << " bfmi=" << posterior.health.min_bfmi
              << " divergences=" << posterior.health.divergences
              << " depth_hits=" << posterior.health.depth_hits << "\n";
  const auto predictive = quadra::sampling::simulate_posterior(
      posterior, {"rep_group_mean[0]", "rep_group_mean[1]",
                  "rep_group_mean[2]"},
      [](const std::vector<double> &q, std::mt19937_64 &rng) {
        std::normal_distribution<double> standard_normal;
        constexpr double tau = 0.4;
        const double sigma = std::exp(q[1]);
        std::vector<double> replicated(3);
        for (std::size_t group = 0; group < replicated.size(); ++group)
          replicated[group] =
              q[0] + tau * q[2 + group] + sigma * standard_normal(rng);
        return replicated;
      },
      simulation_options);

  double posterior_mean_mu = 0.0;
  for (const auto &chain : posterior.fit.chains)
    for (const auto &draw : chain.draws) posterior_mean_mu += draw[0];
  posterior_mean_mu /= posterior.total_draws();

  std::ostringstream draws_csv, parameter_csv, chain_csv, predictive_csv,
      summary_csv;
  quadra::sampling::write_posterior_draws_csv(draws_csv, posterior);
  quadra::sampling::write_parameter_diagnostics_csv(parameter_csv, posterior);
  quadra::sampling::write_chain_diagnostics_csv(chain_csv, posterior);
  quadra::sampling::write_posterior_predictive_csv(predictive_csv, predictive);
  quadra::sampling::write_nuts_summary_csv(summary_csv, posterior);

  if (!posterior.health.passed || posterior.total_draws() != 1600 ||
      predictive.draws.size() != 100 ||
      std::abs(posterior_mean_mu - true_mu) > 0.35 ||
      draws_csv.str().find("mu") == std::string::npos ||
      parameter_csv.str().find("bulk_ess") == std::string::npos ||
      chain_csv.str().find("bfmi") == std::string::npos ||
      predictive_csv.str().find("rep_group_mean") == std::string::npos ||
      summary_csv.str().find("health,PASS") == std::string::npos) {
    std::cerr << "FAIL: end-to-end sampling pipeline\n";
    return 1;
  }
  std::cout << "PASS: simulate -> fit -> sample -> diagnose -> predict -> "
               "serialize\n"
            << "  fitted mu: " << mode.par[0] << "\n"
            << "  posterior mean mu: " << posterior_mean_mu << "\n"
            << "  max R-hat: " << posterior.health.max_rhat << "\n"
            << "  min bulk ESS: " << posterior.health.min_bulk_ess << "\n";
  return 0;
}
