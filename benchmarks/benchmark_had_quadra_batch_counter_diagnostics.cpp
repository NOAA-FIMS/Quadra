#include <chrono>
#include <cmath>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../core/had_quadra.hpp"

DECLARE_ADGRAPH()

namespace {

using Clock = std::chrono::steady_clock;
using had::AReal;
using had::Real;

AReal make_objective(const std::vector<AReal> &x) {
  AReal f(0.0);
  const int n = static_cast<int>(x.size());

  for (int i = 0; i < n; ++i) {
    const auto &xi = x[static_cast<size_t>(i)];
    f = f + exp(0.03 * xi * xi) + sin(0.07 * xi) + 0.01 * xi * xi;
  }

  for (int i = 1; i < n; ++i) {
    const AReal diff =
        x[static_cast<size_t>(i)] - x[static_cast<size_t>(i - 1)];

    f = f + 0.2 * diff * diff +
        0.02 * exp(x[static_cast<size_t>(i)] * x[static_cast<size_t>(i - 1)] *
                   0.01);
  }

  for (int i = 2; i < n; ++i) {
    f = f + 0.003 * x[static_cast<size_t>(i)] * x[static_cast<size_t>(i - 1)] *
                x[static_cast<size_t>(i - 2)];
  }

  return f;
}

struct GraphBundle {
  had::ADGraph graph;
  std::vector<AReal> x;
  AReal f;
};

GraphBundle build_graph(int n) {
  GraphBundle b;
  had::g_ADGraph = &b.graph;

  b.x.reserve(static_cast<size_t>(n));

  for (int i = 0; i < n; ++i) {
    b.x.emplace_back(0.2 + 0.01 * static_cast<double>(i));
  }

  b.f = make_objective(b.x);
  had::g_ADGraph->vertices[b.f.varId].w = 1.0;
  had::PropagateAdjoint();

  return b;
}

void seed_batch(GraphBundle &b, int K) {
  const int n = static_cast<int>(b.x.size());

  had::ResizeDirectionalBatch(K);

  for (int k = 0; k < K; ++k) {
    for (int i = 0; i < n; ++i) {
      const double phase = static_cast<double>((i + 1) * (k + 1));
      const double v = std::sin(0.013 * phase) + 0.5 * std::cos(0.017 * phase);

      had::SetARealDotBatch(b.x[static_cast<size_t>(i)], k, v);
    }
  }
}

template <class Fn> double elapsed_ms(Fn &&fn) {
  const auto start = Clock::now();
  volatile double sink = fn();
  (void)sink;
  const auto end = Clock::now();

  return std::chrono::duration<double, std::milli>(end - start).count();
}

struct Row {
  double total_ms = 0.0;
  std::uint64_t queries = 0;
  std::uint64_t pushdots = 0;
  std::uint64_t inserts = 0;
};

Row run_case(int n, int K) {
  GraphBundle b = build_graph(n);
  had::g_ADGraph = &b.graph;

  seed_batch(b, K);

  Row out;

  out.total_ms = elapsed_ms([&]() {
    had::PropagateAdjointDirectionalBatch();

    // Touch a result so the compiler cannot discard the sweep.
    return had::GetAdjointDotBatch(b.x.back(), b.x.front(), 0);
  });

  out.queries = had::g_batch_query_count;
  out.pushdots = had::g_batch_pushdot_count;
  out.inserts = had::g_batch_insert_count;

  return out;
}

} // namespace

int main(int argc, char **argv) {
  int reps = 1;
  if (argc > 1) {
    reps = std::stoi(argv[1]);
  }

  const std::vector<int> n_values = {100, 250, 500, 1000};
  const std::vector<int> K_values = {1, 2, 4, 8, 16};

  std::cout << "HAD Quadra batch counter diagnostics\n";
  std::cout << "reps per case = " << reps << "\n\n";

  std::cout << std::setw(8) << "n" << std::setw(8) << "K" << std::setw(14)
            << "total ms" << std::setw(14) << "queries" << std::setw(14)
            << "pushdots" << std::setw(14) << "inserts" << std::setw(14)
            << "q/K" << std::setw(14) << "pd/K" << std::setw(14) << "ins/K"
            << "\n";

  std::cout << std::scientific << std::setprecision(6);

  for (int n : n_values) {
    for (int K : K_values) {
      double total_ms = 0.0;
      std::uint64_t queries = 0;
      std::uint64_t pushdots = 0;
      std::uint64_t inserts = 0;

      for (int r = 0; r < reps; ++r) {
        const Row row = run_case(n, K);
        total_ms += row.total_ms;
        queries = row.queries;
        pushdots = row.pushdots;
        inserts = row.inserts;
      }

      total_ms /= static_cast<double>(reps);

      std::cout << std::setw(8) << n << std::setw(8) << K << std::setw(14)
                << total_ms << std::setw(14) << queries << std::setw(14)
                << pushdots << std::setw(14) << inserts << std::setw(14)
                << static_cast<double>(queries) / K << std::setw(14)
                << static_cast<double>(pushdots) / K << std::setw(14)
                << static_cast<double>(inserts) / K << "\n";
    }
  }

  std::cout << "\nBenchmark complete.\n";
  return 0;
}
