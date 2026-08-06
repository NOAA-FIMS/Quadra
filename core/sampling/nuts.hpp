#pragma once

#include "../autodiff.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

namespace quadra {
namespace sampling {

struct LogDensityEvaluation {
  double log_density = -std::numeric_limits<double>::infinity();
  std::vector<double> gradient;
  bool finite = false;
};

template <class LogDensity>
LogDensityEvaluation evaluate_ad_log_density(LogDensity &log_density,
                                             const std::vector<double> &q) {
  had::ADGraph graph;
  had::g_ADGraph = &graph;
  std::vector<AD> q_ad;
  q_ad.reserve(q.size());
  for (double value : q) q_ad.emplace_back(value);

  LogDensityEvaluation out;
  out.gradient.resize(q.size());
  try {
    const AD target = log_density(q_ad);
    out.log_density = value_of(target);
    if (std::isfinite(out.log_density)) {
      had::SetAdjoint(target, 1.0);
      had::PropagateFirstOrderAdjoint();
      out.finite = true;
      for (std::size_t i = 0; i < q.size(); ++i) {
        out.gradient[i] = had::GetAdjoint(q_ad[i]);
        out.finite = out.finite && std::isfinite(out.gradient[i]);
      }
    }
  } catch (const std::exception &) {
    out.finite = false;
  }
  had::g_ADGraph = nullptr;
  return out;
}

struct NutsOptions {
  int warmup = 500;
  int samples = 500;
  int max_tree_depth = 10;
  double target_acceptance = 0.8;
  double initial_step_size = 0.0;
  double divergence_threshold = 1000.0;
  bool adapt_diagonal_mass = true;
  std::uint64_t seed = 5489u;
};

struct NutsDiagnostics {
  int leapfrog_steps = 0;
  int divergences = 0;
  int max_depth_hits = 0;
  double mean_acceptance = 0.0;
  double step_size = 0.0;
  std::vector<double> inverse_mass;
};

struct NutsResult {
  std::vector<std::vector<double>> draws;
  std::vector<double> log_density;
  NutsDiagnostics diagnostics;
};

namespace detail {

inline double dot(const std::vector<double> &a,
                  const std::vector<double> &b) {
  double value = 0.0;
  for (std::size_t i = 0; i < a.size(); ++i) value += a[i] * b[i];
  return value;
}

inline double kinetic(const std::vector<double> &momentum,
                      const std::vector<double> &inverse_mass) {
  double value = 0.0;
  for (std::size_t i = 0; i < momentum.size(); ++i)
    value += momentum[i] * momentum[i] * inverse_mass[i];
  return 0.5 * value;
}

template <class LogDensity>
bool leapfrog(LogDensity &target, std::vector<double> &q,
              std::vector<double> &momentum, LogDensityEvaluation &state,
              const std::vector<double> &inverse_mass, double step) {
  for (std::size_t i = 0; i < q.size(); ++i)
    momentum[i] += 0.5 * step * state.gradient[i];
  for (std::size_t i = 0; i < q.size(); ++i)
    q[i] += step * inverse_mass[i] * momentum[i];
  state = evaluate_ad_log_density(target, q);
  if (!state.finite) return false;
  for (std::size_t i = 0; i < q.size(); ++i)
    momentum[i] += 0.5 * step * state.gradient[i];
  return true;
}

inline bool no_u_turn(const std::vector<double> &q_minus,
                      const std::vector<double> &q_plus,
                      const std::vector<double> &p_minus,
                      const std::vector<double> &p_plus,
                      const std::vector<double> &inverse_mass) {
  std::vector<double> delta(q_minus.size());
  std::vector<double> velocity_minus(q_minus.size());
  std::vector<double> velocity_plus(q_minus.size());
  for (std::size_t i = 0; i < delta.size(); ++i) {
    delta[i] = q_plus[i] - q_minus[i];
    velocity_minus[i] = inverse_mass[i] * p_minus[i];
    velocity_plus[i] = inverse_mass[i] * p_plus[i];
  }
  return dot(delta, velocity_minus) >= 0.0 &&
         dot(delta, velocity_plus) >= 0.0;
}

struct Tree {
  std::vector<double> q_minus, q_plus, p_minus, p_plus, proposal;
  LogDensityEvaluation state_minus, state_plus, proposal_state;
  int valid = 0;
  bool keep_going = true;
  bool divergent = false;
  int leapfrog_steps = 0;
  double acceptance_sum = 0.0;
  int acceptance_count = 0;
};

template <class LogDensity, class Rng>
Tree build_tree(LogDensity &target, const std::vector<double> &q,
                const std::vector<double> &momentum,
                const LogDensityEvaluation &state, double log_slice,
                int direction, int depth, double step_size,
                const std::vector<double> &inverse_mass,
                double initial_joint, double divergence_threshold, Rng &rng) {
  if (depth == 0) {
    Tree tree;
    tree.q_minus = tree.q_plus = tree.proposal = q;
    tree.p_minus = tree.p_plus = momentum;
    tree.state_minus = tree.state_plus = tree.proposal_state = state;
    const bool finite = leapfrog(target, tree.q_minus, tree.p_minus,
                                 tree.state_minus, inverse_mass,
                                 direction * step_size);
    tree.q_plus = tree.proposal = tree.q_minus;
    tree.p_plus = tree.p_minus;
    tree.state_plus = tree.proposal_state = tree.state_minus;
    tree.leapfrog_steps = 1;
    const double joint = finite
                             ? tree.state_minus.log_density -
                                   kinetic(tree.p_minus, inverse_mass)
                             : -std::numeric_limits<double>::infinity();
    const double energy_error = initial_joint - joint;
    tree.valid = finite && log_slice <= joint ? 1 : 0;
    tree.divergent = !finite || !std::isfinite(energy_error) ||
                     energy_error > divergence_threshold;
    tree.keep_going = !tree.divergent &&
                      log_slice < joint + divergence_threshold;
    tree.acceptance_sum =
        finite ? std::min(1.0, std::exp(joint - initial_joint)) : 0.0;
    tree.acceptance_count = 1;
    return tree;
  }

  Tree left = build_tree(target, q, momentum, state, log_slice, direction,
                         depth - 1, step_size, inverse_mass, initial_joint,
                         divergence_threshold, rng);
  if (!left.keep_going) return left;

  Tree right = direction < 0
                   ? build_tree(target, left.q_minus, left.p_minus,
                                left.state_minus, log_slice, direction,
                                depth - 1, step_size, inverse_mass,
                                initial_joint, divergence_threshold, rng)
                   : build_tree(target, left.q_plus, left.p_plus,
                                left.state_plus, log_slice, direction,
                                depth - 1, step_size, inverse_mass,
                                initial_joint, divergence_threshold, rng);
  Tree combined = left;
  if (direction < 0) {
    combined.q_minus = right.q_minus;
    combined.p_minus = right.p_minus;
    combined.state_minus = right.state_minus;
  } else {
    combined.q_plus = right.q_plus;
    combined.p_plus = right.p_plus;
    combined.state_plus = right.state_plus;
  }
  const int total_valid = left.valid + right.valid;
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  if (right.valid > 0 && total_valid > 0 &&
      uniform(rng) < static_cast<double>(right.valid) / total_valid) {
    combined.proposal = right.proposal;
    combined.proposal_state = right.proposal_state;
  }
  combined.valid = total_valid;
  combined.divergent = left.divergent || right.divergent;
  combined.leapfrog_steps += right.leapfrog_steps;
  combined.acceptance_sum += right.acceptance_sum;
  combined.acceptance_count += right.acceptance_count;
  combined.keep_going = left.keep_going && right.keep_going &&
                        no_u_turn(combined.q_minus, combined.q_plus,
                                  combined.p_minus, combined.p_plus,
                                  inverse_mass);
  return combined;
}

template <class LogDensity>
double reasonable_step_size(LogDensity &target, const std::vector<double> &q,
                            const LogDensityEvaluation &state,
                            const std::vector<double> &inverse_mass) {
  std::vector<double> momentum(q.size(), 0.0);
  for (std::size_t i = 0; i < momentum.size(); ++i)
    momentum[i] = 1.0 / std::sqrt(inverse_mass[i]);
  const double initial_joint = state.log_density - kinetic(momentum, inverse_mass);
  double step = 1.0;
  auto acceptance = [&](double candidate) {
    std::vector<double> q_trial = q;
    std::vector<double> p_trial = momentum;
    LogDensityEvaluation trial = state;
    if (!leapfrog(target, q_trial, p_trial, trial, inverse_mass, candidate))
      return 0.0;
    return std::exp(std::min(0.0, trial.log_density -
                                     kinetic(p_trial, inverse_mass) -
                                     initial_joint));
  };
  const bool increase = acceptance(step) > 0.5;
  for (int i = 0; i < 20; ++i) {
    const double next = increase ? step * 2.0 : step * 0.5;
    const double a = acceptance(next);
    step = next;
    if ((increase && a < 0.5) || (!increase && a > 0.5)) break;
  }
  return std::max(1.0e-6, std::min(10.0, step));
}

} // namespace detail

template <class LogDensity>
NutsResult sample_nuts(LogDensity &target, std::vector<double> initial,
                       const NutsOptions &options = {}) {
  if (initial.empty() || options.warmup < 0 || options.samples <= 0 ||
      options.max_tree_depth <= 0 || !(options.target_acceptance > 0.0) ||
      !(options.target_acceptance < 1.0))
    throw std::invalid_argument("sample_nuts: invalid options");

  std::mt19937_64 rng(options.seed);
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  std::normal_distribution<double> normal(0.0, 1.0);
  std::vector<double> inverse_mass(initial.size(), 1.0);
  LogDensityEvaluation current = evaluate_ad_log_density(target, initial);
  if (!current.finite)
    throw std::runtime_error("sample_nuts: initial state is non-finite");

  double step_size = options.initial_step_size > 0.0
                         ? options.initial_step_size
                         : detail::reasonable_step_size(target, initial, current,
                                                        inverse_mass);
  double mu = std::log(10.0 * step_size);
  double log_step_bar = std::log(step_size);
  double h_bar = 0.0;
  std::vector<double> mean(initial.size(), 0.0), m2(initial.size(), 0.0);
  int mass_count = 0;

  NutsResult result;
  result.draws.reserve(static_cast<std::size_t>(options.samples));
  result.log_density.reserve(static_cast<std::size_t>(options.samples));
  const int total_iterations = options.warmup + options.samples;
  double acceptance_total = 0.0;
  int acceptance_iterations = 0;

  for (int iteration = 1; iteration <= total_iterations; ++iteration) {
    std::vector<double> momentum(initial.size());
    for (std::size_t i = 0; i < momentum.size(); ++i)
      momentum[i] = normal(rng) / std::sqrt(inverse_mass[i]);
    const double initial_joint =
        current.log_density - detail::kinetic(momentum, inverse_mass);
    const double log_slice = initial_joint + std::log(uniform(rng));

    detail::Tree tree;
    tree.q_minus = tree.q_plus = tree.proposal = initial;
    tree.p_minus = tree.p_plus = momentum;
    tree.state_minus = tree.state_plus = tree.proposal_state = current;
    tree.valid = 1;
    int depth = 0;
    while (tree.keep_going && depth < options.max_tree_depth) {
      const int direction = uniform(rng) < 0.5 ? -1 : 1;
      detail::Tree extension =
          direction < 0
              ? detail::build_tree(target, tree.q_minus, tree.p_minus,
                                   tree.state_minus, log_slice, direction,
                                   depth, step_size, inverse_mass,
                                   initial_joint,
                                   options.divergence_threshold, rng)
              : detail::build_tree(target, tree.q_plus, tree.p_plus,
                                   tree.state_plus, log_slice, direction,
                                   depth, step_size, inverse_mass,
                                   initial_joint,
                                   options.divergence_threshold, rng);
      if (direction < 0) {
        tree.q_minus = extension.q_minus;
        tree.p_minus = extension.p_minus;
        tree.state_minus = extension.state_minus;
      } else {
        tree.q_plus = extension.q_plus;
        tree.p_plus = extension.p_plus;
        tree.state_plus = extension.state_plus;
      }
      const int total_valid = tree.valid + extension.valid;
      if (extension.valid > 0 && total_valid > 0 &&
          uniform(rng) < static_cast<double>(extension.valid) / total_valid) {
        tree.proposal = extension.proposal;
        tree.proposal_state = extension.proposal_state;
      }
      tree.valid = total_valid;
      tree.divergent = tree.divergent || extension.divergent;
      tree.leapfrog_steps += extension.leapfrog_steps;
      tree.acceptance_sum += extension.acceptance_sum;
      tree.acceptance_count += extension.acceptance_count;
      tree.keep_going = extension.keep_going &&
                        detail::no_u_turn(tree.q_minus, tree.q_plus,
                                          tree.p_minus, tree.p_plus,
                                          inverse_mass);
      ++depth;
    }
    initial = tree.proposal;
    current = tree.proposal_state;
    result.diagnostics.leapfrog_steps += tree.leapfrog_steps;
    result.diagnostics.divergences += tree.divergent ? 1 : 0;
    result.diagnostics.max_depth_hits +=
        depth == options.max_tree_depth ? 1 : 0;
    const double acceptance =
        tree.acceptance_count > 0
            ? tree.acceptance_sum / tree.acceptance_count
            : 0.0;
    acceptance_total += acceptance;
    ++acceptance_iterations;

    if (iteration <= options.warmup) {
      const double eta = 1.0 / (iteration + 10.0);
      h_bar = (1.0 - eta) * h_bar +
              eta * (options.target_acceptance - acceptance);
      const double log_step = mu - std::sqrt(static_cast<double>(iteration)) /
                                       0.05 * h_bar;
      const double weight = std::pow(static_cast<double>(iteration), -0.75);
      log_step_bar = weight * log_step + (1.0 - weight) * log_step_bar;
      step_size = std::exp(log_step);

      if (options.adapt_diagonal_mass) {
        ++mass_count;
        for (std::size_t i = 0; i < initial.size(); ++i) {
          const double delta = initial[i] - mean[i];
          mean[i] += delta / mass_count;
          m2[i] += delta * (initial[i] - mean[i]);
        }
        if (iteration == std::max(20, options.warmup / 2) && mass_count > 1) {
          for (std::size_t i = 0; i < inverse_mass.size(); ++i) {
            const double variance = m2[i] / (mass_count - 1);
            inverse_mass[i] = 1.0 / std::max(1.0e-3, variance);
          }
          step_size = detail::reasonable_step_size(target, initial, current,
                                                    inverse_mass);
          mu = std::log(10.0 * step_size);
          log_step_bar = std::log(step_size);
          h_bar = 0.0;
        }
      }
      if (iteration == options.warmup) step_size = std::exp(log_step_bar);
    } else {
      result.draws.push_back(initial);
      result.log_density.push_back(current.log_density);
    }
  }
  result.diagnostics.mean_acceptance =
      acceptance_total / std::max(1, acceptance_iterations);
  result.diagnostics.step_size = step_size;
  result.diagnostics.inverse_mass = inverse_mass;
  return result;
}

} // namespace sampling
} // namespace quadra
