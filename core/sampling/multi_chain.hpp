#pragma once

#include "nuts.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <thread>
#include <vector>

namespace quadra {
namespace sampling {

struct MultiChainDiagnostics {
  std::vector<double> split_rhat;
  std::vector<double> bulk_ess;
  std::vector<double> tail_ess;
};

struct MultiChainResult {
  std::vector<NutsResult> chains;
  MultiChainDiagnostics diagnostics;
};

namespace detail {

inline double quantile(std::vector<double> values, double probability) {
  if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
  std::sort(values.begin(), values.end());
  const double position = probability * static_cast<double>(values.size() - 1);
  const std::size_t lower = static_cast<std::size_t>(std::floor(position));
  const std::size_t upper = static_cast<std::size_t>(std::ceil(position));
  const double weight = position - lower;
  return (1.0 - weight) * values[lower] + weight * values[upper];
}

inline std::vector<std::vector<double>>
split_parameter_chains(const std::vector<NutsResult> &chains,
                       std::size_t parameter) {
  std::size_t half = std::numeric_limits<std::size_t>::max();
  for (const auto &chain : chains) half = std::min(half, chain.draws.size() / 2);
  if (half < 2) return {};
  std::vector<std::vector<double>> split;
  split.reserve(2 * chains.size());
  for (const auto &chain : chains) {
    std::vector<double> first(half), second(half);
    for (std::size_t i = 0; i < half; ++i) {
      first[i] = chain.draws[i][parameter];
      second[i] = chain.draws[chain.draws.size() - half + i][parameter];
    }
    split.push_back(std::move(first));
    split.push_back(std::move(second));
  }
  return split;
}

inline double sample_variance(const std::vector<double> &x, double mean) {
  double sum = 0.0;
  for (double value : x) {
    const double difference = value - mean;
    sum += difference * difference;
  }
  return x.size() > 1 ? sum / static_cast<double>(x.size() - 1) : 0.0;
}

inline double split_rhat(const std::vector<std::vector<double>> &split) {
  if (split.empty()) return std::numeric_limits<double>::quiet_NaN();
  const std::size_t m = split.size();
  const std::size_t n = split.front().size();
  std::vector<double> means(m), variances(m);
  for (std::size_t chain = 0; chain < m; ++chain) {
    for (double value : split[chain]) means[chain] += value;
    means[chain] /= n;
    variances[chain] = sample_variance(split[chain], means[chain]);
  }
  double mean_of_means = 0.0;
  double within = 0.0;
  for (std::size_t chain = 0; chain < m; ++chain) {
    mean_of_means += means[chain];
    within += variances[chain];
  }
  mean_of_means /= m;
  within /= m;
  double between_sum = 0.0;
  for (double mean : means) {
    const double difference = mean - mean_of_means;
    between_sum += difference * difference;
  }
  const double between =
      static_cast<double>(n) * between_sum / static_cast<double>(m - 1);
  if (!(within > 0.0)) return between == 0.0 ? 1.0 :
      std::numeric_limits<double>::infinity();
  const double variance =
      (static_cast<double>(n - 1) / n) * within + between / n;
  return std::sqrt(variance / within);
}

inline double effective_sample_size(
    const std::vector<std::vector<double>> &split) {
  if (split.empty()) return std::numeric_limits<double>::quiet_NaN();
  const std::size_t m = split.size();
  const std::size_t n = split.front().size();
  std::vector<double> means(m), variances(m);
  for (std::size_t chain = 0; chain < m; ++chain) {
    for (double value : split[chain]) means[chain] += value;
    means[chain] /= n;
    variances[chain] = sample_variance(split[chain], means[chain]);
  }
  double within = 0.0;
  double mean_of_means = 0.0;
  for (std::size_t chain = 0; chain < m; ++chain) {
    within += variances[chain];
    mean_of_means += means[chain];
  }
  within /= m;
  mean_of_means /= m;
  double between_sum = 0.0;
  for (double mean : means) {
    const double difference = mean - mean_of_means;
    between_sum += difference * difference;
  }
  const double between =
      static_cast<double>(n) * between_sum / static_cast<double>(m - 1);
  const double variance =
      (static_cast<double>(n - 1) / n) * within + between / n;
  if (!(variance > 0.0)) return static_cast<double>(m * n);

  std::vector<double> rho(n, 0.0);
  rho[0] = 1.0;
  for (std::size_t lag = 1; lag < n; ++lag) {
    double autocovariance = 0.0;
    for (std::size_t chain = 0; chain < m; ++chain) {
      double chain_covariance = 0.0;
      for (std::size_t i = 0; i + lag < n; ++i)
        chain_covariance +=
            (split[chain][i] - means[chain]) *
            (split[chain][i + lag] - means[chain]);
      autocovariance += chain_covariance / static_cast<double>(n - lag);
    }
    autocovariance /= m;
    rho[lag] = 1.0 - (within - autocovariance) / variance;
  }
  double correlation_sum = 0.0;
  double previous_pair = std::numeric_limits<double>::infinity();
  for (std::size_t lag = 1; lag + 1 < n; lag += 2) {
    double pair = rho[lag] + rho[lag + 1];
    if (pair < 0.0) break;
    pair = std::min(pair, previous_pair);
    correlation_sum += pair;
    previous_pair = pair;
  }
  const double ess = static_cast<double>(m * n) /
                     std::max(1.0, 1.0 + 2.0 * correlation_sum);
  return std::min(static_cast<double>(m * n), ess);
}

inline std::vector<std::vector<double>> indicator_chains(
    const std::vector<std::vector<double>> &split, double threshold,
    bool lower) {
  auto indicators = split;
  for (auto &chain : indicators)
    for (double &value : chain)
      value = lower ? (value <= threshold ? 1.0 : 0.0)
                    : (value >= threshold ? 1.0 : 0.0);
  return indicators;
}

} // namespace detail

inline MultiChainDiagnostics
compute_multi_chain_diagnostics(const std::vector<NutsResult> &chains) {
  if (chains.size() < 2 || chains.front().draws.empty())
    throw std::invalid_argument(
        "multi-chain diagnostics require at least two nonempty chains");
  const std::size_t parameters = chains.front().draws.front().size();
  MultiChainDiagnostics out;
  out.split_rhat.resize(parameters);
  out.bulk_ess.resize(parameters);
  out.tail_ess.resize(parameters);
  for (std::size_t parameter = 0; parameter < parameters; ++parameter) {
    const auto split = detail::split_parameter_chains(chains, parameter);
    if (split.empty())
      throw std::invalid_argument(
          "multi-chain diagnostics require at least four draws per chain");
    std::vector<double> pooled;
    for (const auto &chain : split)
      pooled.insert(pooled.end(), chain.begin(), chain.end());
    const double lower = detail::quantile(pooled, 0.05);
    const double upper = detail::quantile(pooled, 0.95);
    out.split_rhat[parameter] = detail::split_rhat(split);
    out.bulk_ess[parameter] = detail::effective_sample_size(split);
    out.tail_ess[parameter] = std::min(
        detail::effective_sample_size(
            detail::indicator_chains(split, lower, true)),
        detail::effective_sample_size(
            detail::indicator_chains(split, upper, false)));
  }
  return out;
}

template <class TargetFactory>
MultiChainResult sample_nuts_chains(
    TargetFactory target_factory,
    const std::vector<std::vector<double>> &initial_states,
    const NutsOptions &options = {}, bool parallel = true) {
  if (initial_states.size() < 2)
    throw std::invalid_argument("sample_nuts_chains requires at least two chains");
  MultiChainResult out;
  out.chains.resize(initial_states.size());
  std::vector<std::exception_ptr> errors(initial_states.size());
  auto run = [&](std::size_t chain) {
    try {
      auto target = target_factory(chain);
      NutsOptions chain_options = options;
      chain_options.seed +=
          UINT64_C(0x9e3779b97f4a7c15) * static_cast<std::uint64_t>(chain + 1);
      out.chains[chain] =
          sample_nuts(target, initial_states[chain], chain_options);
    } catch (...) {
      errors[chain] = std::current_exception();
    }
  };
  if (parallel) {
    std::vector<std::thread> workers;
    workers.reserve(initial_states.size());
    for (std::size_t chain = 0; chain < initial_states.size(); ++chain)
      workers.emplace_back(run, chain);
    for (auto &worker : workers) worker.join();
  } else {
    for (std::size_t chain = 0; chain < initial_states.size(); ++chain)
      run(chain);
  }
  for (const auto &error : errors)
    if (error) std::rethrow_exception(error);
  out.diagnostics = compute_multi_chain_diagnostics(out.chains);
  return out;
}

} // namespace sampling
} // namespace quadra
