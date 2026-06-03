#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <vector>

#include "../core/had/batch_directional_flat_accumulator.hpp"

namespace {

using Clock = std::chrono::steady_clock;

double ms_between(const Clock::time_point& a, const Clock::time_point& b) {
  return std::chrono::duration<double, std::milli>(b - a).count();
}

struct Row {
  int slots = 0;
  int K = 0;
  int updates = 0;
  double flat_ms = 0.0;
  double map_ms = 0.0;
  double speedup = 0.0;
  double max_diff = 0.0;
};

Row run_case(int slots, int K, int updates, int reps) {
  std::vector<int> update_slot(static_cast<std::size_t>(updates));
  std::vector<int> update_dir(static_cast<std::size_t>(updates));
  std::vector<double> update_value(static_cast<std::size_t>(updates));

  std::mt19937 rng(12345);
  std::uniform_int_distribution<int> slot_dist(0, slots - 1);
  std::uniform_int_distribution<int> dir_dist(0, K - 1);
  std::normal_distribution<double> val_dist(0.0, 1.0);

  for (int i = 0; i < updates; ++i) {
    update_slot[static_cast<std::size_t>(i)] = slot_dist(rng);
    update_dir[static_cast<std::size_t>(i)] = dir_dist(rng);
    update_value[static_cast<std::size_t>(i)] = val_dist(rng);
  }

  had::BatchDirectionalFlatAccumulator flat(
      static_cast<std::size_t>(K),
      static_cast<std::size_t>(slots));

  std::vector<std::map<int, double>> mapped(static_cast<std::size_t>(K));

  const auto f0 = Clock::now();
  for (int r = 0; r < reps; ++r) {
    flat.Clear();

    for (int i = 0; i < updates; ++i) {
      flat.Add(static_cast<std::size_t>(update_dir[static_cast<std::size_t>(i)]),
               static_cast<std::size_t>(update_slot[static_cast<std::size_t>(i)]),
               update_value[static_cast<std::size_t>(i)]);
    }
  }
  const auto f1 = Clock::now();

  const auto m0 = Clock::now();
  for (int r = 0; r < reps; ++r) {
    for (auto& tree : mapped) {
      tree.clear();
    }

    for (int i = 0; i < updates; ++i) {
      mapped[static_cast<std::size_t>(update_dir[static_cast<std::size_t>(i)])]
            [update_slot[static_cast<std::size_t>(i)]] +=
          update_value[static_cast<std::size_t>(i)];
    }
  }
  const auto m1 = Clock::now();

  double max_diff = 0.0;
  for (int k = 0; k < K; ++k) {
    for (int slot = 0; slot < slots; ++slot) {
      const double flat_v =
          flat(static_cast<std::size_t>(k), static_cast<std::size_t>(slot));

      double map_v = 0.0;
      const auto it = mapped[static_cast<std::size_t>(k)].find(slot);
      if (it != mapped[static_cast<std::size_t>(k)].end()) {
        map_v = it->second;
      }

      max_diff = std::max(max_diff, std::abs(flat_v - map_v));
    }
  }

  Row row;
  row.slots = slots;
  row.K = K;
  row.updates = updates;
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

  std::cout << "Batch directional flat accumulator benchmark\n";
  std::cout << "reps per case = " << reps << "\n\n";

  std::cout << std::setw(10) << "slots"
            << std::setw(8) << "K"
            << std::setw(12) << "updates"
            << std::setw(14) << "flat ms"
            << std::setw(14) << "map ms"
            << std::setw(14) << "speedup"
            << std::setw(14) << "max diff"
            << "\n";

  std::cout << std::scientific << std::setprecision(6);

  const std::vector<int> slots_cases = {199, 499, 999};
  const int K = 5;

  for (int slots : slots_cases) {
    const int updates = 12 * slots;
    const Row row = run_case(slots, K, updates, reps);

    std::cout << std::setw(10) << row.slots
              << std::setw(8) << row.K
              << std::setw(12) << row.updates
              << std::setw(14) << row.flat_ms
              << std::setw(14) << row.map_ms
              << std::setw(14) << row.speedup
              << std::setw(14) << row.max_diff
              << "\n";
  }

  std::cout << "\nBenchmark complete.\n";
  return 0;
}
