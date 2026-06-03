#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <vector>

#include "../core/had/batch_directional_flat_accumulator.hpp"
#include "../core/had/intermediate_edge_slot_registry.hpp"

namespace {

using Clock = std::chrono::steady_clock;

double ms_between(const Clock::time_point& a, const Clock::time_point& b) {
  return std::chrono::duration<double, std::milli>(b - a).count();
}

struct Update {
  std::uint32_t a = 0;
  std::uint32_t b = 0;
  int direction = 0;
  double value = 0.0;
};

struct Row {
  int vertices = 0;
  int edges = 0;
  int updates = 0;
  int K = 0;
  double registry_build_ms = 0.0;
  double flat_ms = 0.0;
  double map_ms = 0.0;
  double speedup = 0.0;
  double max_diff = 0.0;
};

std::uint64_t pack(std::uint32_t a, std::uint32_t b) {
  const std::uint32_t outer = std::max(a, b);
  const std::uint32_t inner = std::min(a, b);
  return (static_cast<std::uint64_t>(outer) << 32) |
         static_cast<std::uint64_t>(inner);
}

Row run_case(int vertices, int target_edges, int K, int reps) {
  std::mt19937 rng(12345);
  std::uniform_int_distribution<int> vertex_dist(0, vertices - 1);
  std::uniform_int_distribution<int> direction_dist(0, K - 1);
  std::normal_distribution<double> value_dist(0.0, 1.0);

  std::vector<Update> updates;
  updates.reserve(static_cast<std::size_t>(target_edges * 4));

  had::IntermediateEdgeSlotRegistry registry;

  const auto rb0 = Clock::now();

  while (static_cast<int>(registry.size()) < target_edges) {
    const auto a = static_cast<std::uint32_t>(vertex_dist(rng));
    const auto b = static_cast<std::uint32_t>(vertex_dist(rng));
    if (a == b) continue;

    registry.GetOrCreate(a, b);
  }

  const auto rb1 = Clock::now();

  const auto& edges = registry.edges();

  for (const auto& edge : edges) {
    for (int pass = 0; pass < 4; ++pass) {
      updates.push_back(Update{
          edge.outer,
          edge.inner,
          direction_dist(rng),
          value_dist(rng)});
    }
  }

  had::BatchDirectionalFlatAccumulator flat(
      static_cast<std::size_t>(K),
      registry.size());

  std::vector<std::map<std::uint64_t, double>> mapped(
      static_cast<std::size_t>(K));

  const auto f0 = Clock::now();

  for (int r = 0; r < reps; ++r) {
    flat.Clear();

    for (const auto& u : updates) {
      std::size_t slot = 0;
      if (!registry.TryGet(u.a, u.b, slot)) {
        throw std::runtime_error("missing registry slot");
      }

      flat.Add(static_cast<std::size_t>(u.direction), slot, u.value);
    }
  }

  const auto f1 = Clock::now();

  const auto m0 = Clock::now();

  for (int r = 0; r < reps; ++r) {
    for (auto& direction_map : mapped) {
      direction_map.clear();
    }

    for (const auto& u : updates) {
      mapped[static_cast<std::size_t>(u.direction)][pack(u.a, u.b)] += u.value;
    }
  }

  const auto m1 = Clock::now();

  double max_diff = 0.0;

  for (const auto& edge : edges) {
    for (int k = 0; k < K; ++k) {
      const double flat_value =
          flat(static_cast<std::size_t>(k), edge.slot);

      double map_value = 0.0;
      const auto key = pack(edge.outer, edge.inner);
      const auto found = mapped[static_cast<std::size_t>(k)].find(key);
      if (found != mapped[static_cast<std::size_t>(k)].end()) {
        map_value = found->second;
      }

      max_diff = std::max(max_diff, std::abs(flat_value - map_value));
    }
  }

  Row row;
  row.vertices = vertices;
  row.edges = static_cast<int>(registry.size());
  row.updates = static_cast<int>(updates.size());
  row.K = K;
  row.registry_build_ms = ms_between(rb0, rb1);
  row.flat_ms = ms_between(f0, f1) / static_cast<double>(reps);
  row.map_ms = ms_between(m0, m1) / static_cast<double>(reps);
  row.speedup = row.map_ms / row.flat_ms;
  row.max_diff = max_diff;
  return row;
}

}  // namespace

int main(int argc, char** argv) {
  int reps = 100;
  if (argc > 1) {
    reps = std::stoi(argv[1]);
  }

  const int K = 5;

  std::cout << "Intermediate edge slot registry benchmark\n";
  std::cout << "reps per case = " << reps << "\n";
  std::cout << "K = " << K << "\n\n";

  std::cout << std::setw(10) << "vertices"
            << std::setw(10) << "edges"
            << std::setw(12) << "updates"
            << std::setw(14) << "reg ms"
            << std::setw(14) << "flat ms"
            << std::setw(14) << "map ms"
            << std::setw(14) << "speedup"
            << std::setw(14) << "max diff"
            << "\n";

  std::cout << std::scientific << std::setprecision(6);

  for (const auto& cfg : {std::pair<int, int>{2606, 16000},
                          std::pair<int, int>{6506, 40000},
                          std::pair<int, int>{13006, 80000}}) {
    const Row row = run_case(cfg.first, cfg.second, K, reps);

    std::cout << std::setw(10) << row.vertices
              << std::setw(10) << row.edges
              << std::setw(12) << row.updates
              << std::setw(14) << row.registry_build_ms
              << std::setw(14) << row.flat_ms
              << std::setw(14) << row.map_ms
              << std::setw(14) << row.speedup
              << std::setw(14) << row.max_diff
              << "\n";
  }

  std::cout << "\nBenchmark complete.\n";
  return 0;
}
