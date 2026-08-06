#include "../include/quadra/sampling.hpp"

#include <cmath>
#include <iostream>
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
  return 0;
}
