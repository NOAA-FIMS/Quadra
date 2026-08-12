#include <Eigen/Dense>

#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../core/laplace/exact_gradient_workspace.hpp"

DECLARE_ADGRAPH()

namespace {

using Clock = std::chrono::steady_clock;

had::AReal objective(const std::vector<had::AReal> &theta,
                     const std::vector<had::AReal> &u) {
  had::AReal out(0.0);
  for (std::size_t i = 0; i < u.size(); ++i) {
    const had::AReal centered = u[i] - theta[i % theta.size()];
    out = out + 0.5 * centered * centered + exp(theta[0]) * u[i] * u[i];
  }
  return out;
}

} // namespace

int main(int argc, char **argv) {
  const int random_dim = argc > 1 ? std::atoi(argv[1]) : 1000;
  const int directions = argc > 2 ? std::atoi(argv[2]) : 32;
  const std::size_t budget_mb =
      argc > 3 ? static_cast<std::size_t>(std::strtoull(argv[3], nullptr, 10))
               : 64;
  if (random_dim <= 0 || directions <= 0 || budget_mb == 0) {
    std::cerr << "usage: benchmark_adaptive_directional_batch "
                 "[random_dim] [directions] [budget_mb]\n";
    return 2;
  }

  std::vector<had::AReal> theta(static_cast<std::size_t>(directions));
  std::vector<had::AReal> u(static_cast<std::size_t>(random_dim));
  quadra::laplace::ExactGradientWorkspace workspace;
  workspace.Build(
      [&]() {
        for (int j = 0; j < directions; ++j)
          theta[static_cast<std::size_t>(j)] = had::AReal(0.01 * (j + 1));
        for (int i = 0; i < random_dim; ++i)
          u[static_cast<std::size_t>(i)] = had::AReal(0.001 * (i + 1));
        return objective(theta, u);
      },
      &theta, &u);
  workspace.PropagateBaseAdjoint();

  std::vector<quadra::laplace::SparseHdotPatternEntry> pattern;
  pattern.reserve(static_cast<std::size_t>(random_dim));
  for (int i = 0; i < random_dim; ++i)
    pattern.emplace_back(i, i);

  quadra::laplace::AdaptiveDirectionalBatchOptions options;
  options.memory_budget_bytes = budget_mb * 1024u * 1024u;
  options.maximum_batch_size = 64;

  const auto start = Clock::now();
  const auto result = workspace.TraceTermsSelectedInverseAdaptive(
      static_cast<std::size_t>(directions),
      [&](std::size_t k, Eigen::VectorXd &td, Eigen::VectorXd &ud) {
        td = Eigen::VectorXd::Zero(directions);
        td[static_cast<Eigen::Index>(k)] = 1.0;
        ud = Eigen::VectorXd::Zero(random_dim);
      },
      [](int row, int col) { return row == col ? 1.0 : 0.0; }, pattern,
      options);
  const auto end = Clock::now();
  const double elapsed_ms =
      std::chrono::duration<double, std::milli>(end - start).count();

  std::cout << "vertices,budget_mb,batch_size,batches,estimated_lane_bytes,"
               "peak_tracked_bytes,runtime_ms,trace_sum\n";
  std::cout << workspace.HadWorkspace().VertexCount() << ',' << budget_mb << ','
            << result.plan.batch_size << ',' << result.batches_executed << ','
            << result.plan.estimated_bytes_per_direction << ','
            << result.peak_tracked_graph_bytes << ',' << std::fixed
            << std::setprecision(3) << elapsed_ms << ','
            << result.trace_terms.sum() << '\n';
}
