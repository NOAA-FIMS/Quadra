#include "../core/laplace/adaptive_directional_batch.hpp"

#include <cassert>
#include <iostream>
#include <limits>

DECLARE_ADGRAPH()

int main() {
  had::ADGraph graph;
  graph.vertices.reserve(100);
  for (std::size_t i = 0; i < 100; ++i) {
    graph.vertices.emplace_back(static_cast<had::VertexId>(i));
  }

  const std::size_t base =
      had::MeasureADGraphMemory(graph).total_tracked_reserved_bytes;
  const std::size_t lane =
      quadra::laplace::EstimateDirectionalLaneBytes(graph.vertices.size());

  quadra::laplace::AdaptiveDirectionalBatchOptions options;
  options.memory_budget_bytes = base + 4 * lane;
  options.maximum_batch_size = 64;
  auto plan =
      quadra::laplace::PlanAdaptiveDirectionalBatch(graph, 20, options);
  assert(plan.batch_size == 4);
  assert(plan.constrained_by_memory);

  options.memory_budget_bytes = base + 100 * lane;
  plan = quadra::laplace::PlanAdaptiveDirectionalBatch(graph, 20, options);
  assert(plan.batch_size == 20);
  assert(!plan.constrained_by_memory);

  options.memory_budget_bytes = 0;
  plan = quadra::laplace::PlanAdaptiveDirectionalBatch(graph, 20, options);
  assert(plan.batch_size == 1);
  assert(plan.constrained_by_memory);

  plan = quadra::laplace::PlanAdaptiveDirectionalBatch(graph, 0, options);
  assert(plan.batch_size == 0);

  std::cout << "adaptive directional batch policy tests passed\n";
}
