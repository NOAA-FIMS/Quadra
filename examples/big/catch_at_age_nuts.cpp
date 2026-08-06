#include "../../include/quadra/sampling.hpp"
#include "catch_at_age_shared.hpp"

#include <cmath>
#include <iostream>
#include <numeric>
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
  options.warmup = 100;
  options.samples = 100;
  options.max_tree_depth = 8;
  options.target_acceptance = 0.85;
  options.seed = 20260806;
  const auto fit = quadra::sampling::sample_nuts(posterior, initial, options);

  double mean_log_m = 0.0;
  for (const auto &draw : fit.draws) mean_log_m += draw[1];
  mean_log_m /= fit.draws.size();
  std::cout << "Quadra non-centered joint AD-NUTS catch-at-age example\n"
            << "draws: " << fit.draws.size() << "\n"
            << "mean log_M: " << mean_log_m << "\n"
            << "mean acceptance: " << fit.diagnostics.mean_acceptance << "\n"
            << "step size: " << fit.diagnostics.step_size << "\n"
            << "leapfrog steps: " << fit.diagnostics.leapfrog_steps << "\n"
            << "divergences: " << fit.diagnostics.divergences << "\n"
            << "max-depth hits: " << fit.diagnostics.max_depth_hits << "\n";
  return fit.diagnostics.divergences > 10 ||
                 fit.diagnostics.max_depth_hits > options.samples
             ? 1
             : 0;
}
