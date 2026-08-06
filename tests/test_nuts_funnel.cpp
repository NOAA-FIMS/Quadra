#include "../include/quadra/sampling.hpp"

#include <cmath>
#include <iostream>
#include <vector>

DECLARE_ADGRAPH();

struct CenteredFunnel {
  template <class T> T operator()(const std::vector<T> &q) const {
    using std::exp;
    const T log_scale = q[0];
    T log_density = -T(0.5 / 9.0) * log_scale * log_scale;
    for (std::size_t i = 1; i < q.size(); ++i)
      log_density -=
          T(0.5) * q[i] * q[i] * exp(-log_scale) + T(0.5) * log_scale;
    return log_density;
  }
};

struct NoncenteredFunnel {
  template <class T> T operator()(const std::vector<T> &q) const {
    T log_density = -T(0.5 / 9.0) * q[0] * q[0];
    for (std::size_t i = 1; i < q.size(); ++i)
      log_density -= T(0.5) * q[i] * q[i];
    return log_density;
  }
};

int main() {
  quadra::sampling::NutsOptions options;
  options.warmup = 500;
  options.samples = 500;
  options.max_tree_depth = 9;
  options.target_acceptance = 0.85;
  options.adapt_dense_mass = true;
  options.seed = 20260810;

  std::vector<std::vector<double>> initial_states(4,
                                                  std::vector<double>(6, 0.0));
  initial_states[0][0] = -8.0;
  initial_states[1][0] = -3.0;
  initial_states[2][0] = 3.0;
  initial_states[3][0] = 8.0;

  const auto centered = quadra::sampling::sample_nuts_chains(
      [](std::size_t) { return CenteredFunnel{}; }, initial_states, options,
      true);
  const auto noncentered = quadra::sampling::sample_nuts_chains(
      [](std::size_t) { return NoncenteredFunnel{}; }, initial_states, options,
      true);
  const auto centered_health = quadra::sampling::assess_nuts_health(centered);
  const auto noncentered_health =
      quadra::sampling::assess_nuts_health(noncentered);

  std::cout << "centered funnel: "
            << (centered_health.passed ? "PASS" : "FLAGGED")
            << " rhat=" << centered_health.max_rhat
            << " bulk_ess=" << centered_health.min_bulk_ess
            << " divergences=" << centered_health.divergences << "\n"
            << "noncentered funnel: "
            << (noncentered_health.passed ? "PASS" : "FLAGGED")
            << " rhat=" << noncentered_health.max_rhat
            << " bulk_ess=" << noncentered_health.min_bulk_ess
            << " divergences=" << noncentered_health.divergences << "\n";

  if (centered_health.passed || !noncentered_health.passed) {
    std::cerr << "FAIL: funnel geometry was not diagnosed correctly\n";
    return 1;
  }
  return 0;
}
