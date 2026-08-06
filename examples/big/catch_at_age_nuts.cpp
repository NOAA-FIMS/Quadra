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

struct CatchAtAgePosteriorPredictive {
  const example::CatchAtAgeLaplaceModel *model;

  std::vector<std::string> quantity_names() const {
    std::vector<std::string> names;
    for (int year = 0; year < model->data.n_years; ++year) {
      names.push_back("pred_index[" + std::to_string(year) + "]");
      names.push_back("rep_index[" + std::to_string(year) + "]");
      names.push_back("pred_catch[" + std::to_string(year) + "]");
      names.push_back("rep_catch[" + std::to_string(year) + "]");
      for (int age = 0; age < model->data.n_ages; ++age)
        names.push_back("pred_comp[" + std::to_string(year) + "," +
                        std::to_string(age) + "]");
      for (int age = 0; age < model->data.n_ages; ++age)
        names.push_back("rep_comp_count[" + std::to_string(year) + "," +
                        std::to_string(age) + "]");
    }
    return names;
  }

  std::vector<double> operator()(const std::vector<double> &q,
                                 std::mt19937_64 &rng) const {
    const int years = model->data.n_years;
    const int ages = model->data.n_ages;
    const double R0 = std::exp(q[0]);
    const double M = std::exp(q[1]);
    const double survey_q = std::exp(q[2]);
    const double Fbar = std::exp(q[3]);
    const double sel50 = example::bounded_sel50(q[4], ages);
    const double sel_slope = std::exp(q[5]);
    const double sigma_index = example::positive_floor(q[6], 0.03);
    const double sigma_catch = example::positive_floor(q[7], 0.08);
    const double sigma_rec = example::positive_floor(q[8], 0.05);
    const double concentration = example::positive_floor(q[9], 1.0);
    std::normal_distribution<double> normal;

    std::vector<double> abundance(static_cast<std::size_t>(ages));
    std::vector<double> next_abundance(static_cast<std::size_t>(ages));
    std::vector<double> selectivity(static_cast<std::size_t>(ages));
    std::vector<double> catch_at_age(static_cast<std::size_t>(ages));
    for (int age = 0; age < ages; ++age) {
      selectivity[static_cast<std::size_t>(age)] =
          example::normalized_logistic_selectivity(
              static_cast<double>(age + 1), sel50, sel_slope, ages);
      abundance[static_cast<std::size_t>(age)] = R0 * std::exp(-M * age);
    }

    std::vector<double> out;
    out.reserve(static_cast<std::size_t>(years * (4 + 2 * ages)));
    for (int year = 0; year < years; ++year) {
      double vulnerable_biomass = 0.0;
      double total_catch = 0.0;
      for (int age = 0; age < ages; ++age) {
        const double age_value = static_cast<double>(age + 1);
        const double weight = age_value * age_value * age_value / 100.0;
        const double fishing_mortality =
            Fbar * selectivity[static_cast<std::size_t>(age)];
        const double total_mortality = M + fishing_mortality;
        const double catch_number =
            abundance[static_cast<std::size_t>(age)] *
            fishing_mortality / total_mortality *
            (1.0 - std::exp(-total_mortality));
        catch_at_age[static_cast<std::size_t>(age)] = catch_number;
        total_catch += catch_number * weight;
        vulnerable_biomass += abundance[static_cast<std::size_t>(age)] *
                              weight *
                              selectivity[static_cast<std::size_t>(age)];
      }
      const double predicted_index = survey_q * vulnerable_biomass + 1.0e-9;
      const double predicted_catch = total_catch + 1.0e-9;
      out.push_back(predicted_index);
      out.push_back(predicted_index * std::exp(sigma_index * normal(rng)));
      out.push_back(predicted_catch);
      out.push_back(predicted_catch * std::exp(sigma_catch * normal(rng)));

      const double catch_sum =
          std::accumulate(catch_at_age.begin(), catch_at_age.end(), 1.0e-12);
      std::vector<double> composition(static_cast<std::size_t>(ages));
      std::vector<double> dirichlet_weights(static_cast<std::size_t>(ages));
      double dirichlet_sum = 0.0;
      for (int age = 0; age < ages; ++age) {
        composition[static_cast<std::size_t>(age)] =
            catch_at_age[static_cast<std::size_t>(age)] / catch_sum;
        out.push_back(composition[static_cast<std::size_t>(age)]);
        std::gamma_distribution<double> gamma(
            std::max(1.0e-12,
                     concentration * composition[static_cast<std::size_t>(age)]),
            1.0);
        dirichlet_weights[static_cast<std::size_t>(age)] = gamma(rng);
        dirichlet_sum += dirichlet_weights[static_cast<std::size_t>(age)];
      }
      if (dirichlet_sum > 0.0 && std::isfinite(dirichlet_sum)) {
        for (double &weight : dirichlet_weights) weight /= dirichlet_sum;
      } else {
        dirichlet_weights = composition;
      }
      std::discrete_distribution<int> age_draw(dirichlet_weights.begin(),
                                               dirichlet_weights.end());
      std::vector<int> replicated_counts(static_cast<std::size_t>(ages), 0);
      const int sample_size = std::accumulate(
          model->data.age_comp_counts[static_cast<std::size_t>(year)].begin(),
          model->data.age_comp_counts[static_cast<std::size_t>(year)].end(), 0);
      for (int fish = 0; fish < sample_size; ++fish)
        ++replicated_counts[static_cast<std::size_t>(age_draw(rng))];
      for (int count : replicated_counts) out.push_back(count);

      std::fill(next_abundance.begin(), next_abundance.end(), 0.0);
      const double recruitment_deviation =
          sigma_rec * q[static_cast<std::size_t>(10 + year)];
      next_abundance[0] =
          R0 * std::exp(recruitment_deviation - 0.5 * sigma_rec * sigma_rec);
      for (int age = 1; age < ages; ++age) {
        const double mortality =
            M + Fbar * selectivity[static_cast<std::size_t>(age - 1)];
        next_abundance[static_cast<std::size_t>(age)] =
            abundance[static_cast<std::size_t>(age - 1)] *
            std::exp(-mortality);
      }
      const double plus_mortality =
          M + Fbar * selectivity[static_cast<std::size_t>(ages - 1)];
      next_abundance[static_cast<std::size_t>(ages - 1)] +=
          abundance[static_cast<std::size_t>(ages - 1)] *
          std::exp(-plus_mortality);
      abundance.swap(next_abundance);
    }
    return out;
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
  options.target_acceptance = 0.85;
  options.adapt_dense_mass = true;
  options.seed = 20260806;
  std::vector<std::string> parameter_names{
      "log_R0",          "log_M",        "log_q",       "log_F",
      "logit_sel50",     "log_sel_slope", "log_sigma_I", "log_sigma_C",
      "log_sigma_rec",   "log_comp_concentration"};
  for (int year = 0; year < model.data.n_years; ++year)
    parameter_names.push_back("z_rec[" + std::to_string(year) + "]");
  quadra::sampling::NutsWorkflowOptions workflow_options;
  workflow_options.sampler = options;
  workflow_options.chains = 4;
  workflow_options.initial_jitter = 0.01;
  workflow_options.initialization_seed = 20260806;
  workflow_options.health.max_rhat = 1.05;
  const auto workflow = quadra::sampling::run_nuts_workflow(
      [&model](std::size_t) {
        return NoncenteredCatchAtAgePosterior{&model};
      },
      initial, parameter_names, workflow_options);
  const auto &fit = workflow.fit;
  CatchAtAgePosteriorPredictive predictive_simulator{&model};
  quadra::sampling::PosteriorSimulationOptions simulation_options;
  simulation_options.thin = 10;
  simulation_options.max_draws = 100;
  simulation_options.seed = 20260813;
  const auto simulations = quadra::sampling::simulate_posterior(
      workflow, predictive_simulator.quantity_names(), predictive_simulator,
      simulation_options);

  int predictive_coverage = 0;
  const std::size_t quantities_per_year =
      static_cast<std::size_t>(4 + 2 * model.data.n_ages);
  for (const auto &draw : simulations.draws) {
    for (int year = 0; year < model.data.n_years; ++year) {
      const std::size_t offset =
          static_cast<std::size_t>(year) * quantities_per_year + 4 +
          static_cast<std::size_t>(model.data.n_ages);
      const double replicated_total = std::accumulate(
          draw.values.begin() + static_cast<std::ptrdiff_t>(offset),
          draw.values.begin() + static_cast<std::ptrdiff_t>(
                                    offset + model.data.n_ages),
          0.0);
      const int observed_total = std::accumulate(
          model.data.age_comp_counts[static_cast<std::size_t>(year)].begin(),
          model.data.age_comp_counts[static_cast<std::size_t>(year)].end(), 0);
      if (replicated_total != observed_total)
        throw std::runtime_error(
            "posterior predictive composition sample size changed");
    }
  }
  auto interval_contains = [&](std::size_t quantity, double observed) {
    std::vector<double> values;
    values.reserve(simulations.draws.size());
    for (const auto &draw : simulations.draws)
      values.push_back(draw.values[quantity]);
    std::sort(values.begin(), values.end());
    const std::size_t lower = static_cast<std::size_t>(0.05 * values.size());
    const std::size_t upper =
        std::min(values.size() - 1,
                 static_cast<std::size_t>(0.95 * values.size()));
    return observed >= values[lower] && observed <= values[upper];
  };
  for (int year = 0; year < model.data.n_years; ++year) {
    const std::size_t offset = static_cast<std::size_t>(year) *
                               quantities_per_year;
    predictive_coverage +=
        interval_contains(offset + 1, model.data.index_obs[year]) ? 1 : 0;
    predictive_coverage +=
        interval_contains(offset + 3, model.data.catch_obs[year]) ? 1 : 0;
  }

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
  const auto &health = workflow.health;
  const std::size_t worst_rhat_parameter =
      static_cast<std::size_t>(worst_rhat - fit.diagnostics.split_rhat.begin());
  std::cout << "Quadra non-centered joint AD-NUTS catch-at-age example\n"
            << "chains: " << fit.chains.size() << "\n"
            << "draws: " << total_draws << "\n"
            << "posterior predictive draws: " << simulations.draws.size()
            << "\n"
            << "90% predictive coverage: " << predictive_coverage << "/"
            << 2 * model.data.n_years << "\n"
            << "mean log_M: " << mean_log_m << "\n"
            << "log_M split R-hat: " << fit.diagnostics.split_rhat[1] << "\n"
            << "log_M bulk ESS: " << fit.diagnostics.bulk_ess[1] << "\n"
            << "log_M tail ESS: " << fit.diagnostics.tail_ess[1] << "\n"
            << "maximum split R-hat: " << *worst_rhat << " ("
            << workflow.parameter_names[worst_rhat_parameter] << ")\n"
            << "minimum bulk ESS: " << *worst_bulk << "\n"
            << "minimum tail ESS: " << *worst_tail << "\n"
            << "minimum BFMI: " << health.min_bfmi << "\n"
            << "divergences: " << divergences << "\n"
            << "max-depth hits: " << depth_hits << "\n"
            << "sampler health: " << (health.passed ? "PASS" : "FAIL")
            << "\n";
  for (std::size_t chain = 0; chain < fit.chains.size(); ++chain) {
    const auto &diagnostics = fit.chains[chain].diagnostics;
    std::cout << "chain " << chain + 1
              << ": acceptance=" << diagnostics.mean_acceptance
              << " step_size=" << diagnostics.step_size
              << " divergences=" << diagnostics.divergences
              << " depth_hits=" << diagnostics.max_depth_hits
              << " BFMI=" << diagnostics.energy_bfmi << "\n";
  }
  return health.passed ? 0 : 1;
}
