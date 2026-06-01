#!/usr/bin/env bash
set -euo pipefail

mkdir -p benchmarks .quadra_patch_backups

if [[ ! -f core/had_quadra.hpp ]]; then
  echo "ERROR: missing core/had_quadra.hpp"
  exit 1
fi

if ! grep -q "PropagateAdjointDirectionalBatch" core/had_quadra.hpp; then
  echo "ERROR: PropagateAdjointDirectionalBatch not found."
  exit 1
fi

cat > benchmarks/benchmark_had_quadra_batch_hotspots.cpp <<'EOF'
#include <cmath>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../core/had_quadra.hpp"

DECLARE_ADGRAPH()

namespace {

using Clock = std::chrono::steady_clock;
using had::AReal;
using had::Real;

AReal make_objective(const std::vector<AReal>& x) {
    AReal f(0.0);
    const int n = static_cast<int>(x.size());

    for (int i = 0; i < n; ++i) {
        const auto& xi = x[static_cast<size_t>(i)];
        f = f + exp(0.03 * xi * xi) + sin(0.07 * xi) + 0.01 * xi * xi;
    }

    for (int i = 1; i < n; ++i) {
        const AReal diff = x[static_cast<size_t>(i)] - x[static_cast<size_t>(i - 1)];
        f = f + 0.2 * diff * diff;
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

void seed_scalar(GraphBundle& b, int k) {
    const int n = static_cast<int>(b.x.size());

    for (int i = 0; i < n; ++i) {
        const double phase = static_cast<double>((i + 1) * (k + 1));
        const double v = std::sin(0.013 * phase) + 0.5 * std::cos(0.017 * phase);

        b.x[static_cast<size_t>(i)].dot = v;
        had::g_ADGraph->vertices[b.x[static_cast<size_t>(i)].varId].dot = v;
    }
}

void seed_batch(GraphBundle& b, int K) {
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

double checksum_scalar(GraphBundle& b) {
    double sum = 0.0;
    const int n = static_cast<int>(b.x.size());
    const int stride = std::max(1, n / 16);

    for (int i = 0; i < n; i += stride) {
        for (int j = 0; j <= i; j += stride) {
            sum += had::GetAdjointDot(b.x[static_cast<size_t>(i)],
                                      b.x[static_cast<size_t>(j)]);
        }
    }

    return sum;
}

double checksum_batch(GraphBundle& b, int K) {
    double sum = 0.0;
    const int n = static_cast<int>(b.x.size());
    const int stride = std::max(1, n / 16);

    for (int k = 0; k < K; ++k) {
        for (int i = 0; i < n; i += stride) {
            for (int j = 0; j <= i; j += stride) {
                sum += had::GetAdjointDotBatch(b.x[static_cast<size_t>(i)],
                                               b.x[static_cast<size_t>(j)],
                                               k);
            }
        }
    }

    return sum;
}

template <class Fn>
double mean_ms(Fn&& fn, int reps) {
    const auto start = Clock::now();

    for (int r = 0; r < reps; ++r) {
        volatile double sink = fn();
        (void)sink;
    }

    const auto end = Clock::now();
    std::chrono::duration<double, std::milli> elapsed = end - start;
    return elapsed.count() / static_cast<double>(reps);
}

struct Row {
    double build_adj_ms = 0.0;
    double scalar_prop_ms = 0.0;
    double scalar_query_ms = 0.0;
    double batch_prop_ms = 0.0;
    double batch_query_ms = 0.0;
    double prop_speedup = 0.0;
};

Row run_case(int n, int K, int reps) {
    Row out;

    out.build_adj_ms = mean_ms([&]() {
        GraphBundle b = build_graph(n);
        return static_cast<double>(b.graph.vertices.size());
    }, reps);

    out.scalar_prop_ms = mean_ms([&]() {
        GraphBundle b = build_graph(n);
        had::g_ADGraph = &b.graph;

        double acc = 0.0;
        for (int k = 0; k < K; ++k) {
            seed_scalar(b, k);
            had::PropagateAdjointDirectional();
            acc += had::GetAdjointDot(b.x.back(), b.x.front());
        }
        return acc;
    }, reps);

    {
        GraphBundle b = build_graph(n);
        had::g_ADGraph = &b.graph;

        for (int k = 0; k < K; ++k) {
            seed_scalar(b, k);
            had::PropagateAdjointDirectional();
        }

        out.scalar_query_ms = mean_ms([&]() {
            return checksum_scalar(b);
        }, reps);
    }

    out.batch_prop_ms = mean_ms([&]() {
        GraphBundle b = build_graph(n);
        had::g_ADGraph = &b.graph;

        seed_batch(b, K);
        had::PropagateAdjointDirectionalBatch();
        return had::GetAdjointDotBatch(b.x.back(), b.x.front(), 0);
    }, reps);

    {
        GraphBundle b = build_graph(n);
        had::g_ADGraph = &b.graph;

        seed_batch(b, K);
        had::PropagateAdjointDirectionalBatch();

        out.batch_query_ms = mean_ms([&]() {
            return checksum_batch(b, K);
        }, reps);
    }

    out.prop_speedup = out.scalar_prop_ms / out.batch_prop_ms;
    return out;
}

}  // namespace

int main(int argc, char** argv) {
    int reps = 10;
    if (argc > 1) {
        reps = std::stoi(argv[1]);
    }

    const std::vector<int> n_values = {100, 250, 500};
    const std::vector<int> K_values = {1, 4, 16};

    std::cout << "HAD Quadra batch propagation hotspot benchmark\n";
    std::cout << "reps per case = " << reps << "\n\n";

    std::cout << std::setw(8) << "n"
              << std::setw(8) << "K"
              << std::setw(14) << "build"
              << std::setw(14) << "s_prop"
              << std::setw(14) << "s_query"
              << std::setw(14) << "b_prop"
              << std::setw(14) << "b_query"
              << std::setw(14) << "prop spd"
              << "\n";

    std::cout << std::scientific << std::setprecision(6);

    for (int n : n_values) {
        for (int K : K_values) {
            const Row r = run_case(n, K, reps);

            std::cout << std::setw(8) << n
                      << std::setw(8) << K
                      << std::setw(14) << r.build_adj_ms
                      << std::setw(14) << r.scalar_prop_ms
                      << std::setw(14) << r.scalar_query_ms
                      << std::setw(14) << r.batch_prop_ms
                      << std::setw(14) << r.batch_query_ms
                      << std::setw(14) << r.prop_speedup
                      << "\n";
        }
    }

    std::cout << "\nBenchmark complete.\n";
    return 0;
}
EOF

cat > run_had_quadra_batch_hotspot_benchmark.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O3 -DNDEBUG}"
REPS="${1:-10}"

EIGEN_INCLUDE=""
if [[ -d external/Eigen ]]; then
  EIGEN_INCLUDE="-Iexternal/Eigen"
elif [[ -d core/eigen ]]; then
  EIGEN_INCLUDE="-Icore/eigen"
fi

mkdir -p build/benchmarks

set -x
"${CXX}" ${CXXFLAGS} ${EIGEN_INCLUDE} -I.   benchmarks/benchmark_had_quadra_batch_hotspots.cpp   -o build/benchmarks/benchmark_had_quadra_batch_hotspots

./build/benchmarks/benchmark_had_quadra_batch_hotspots "${REPS}"
EOF

chmod +x run_had_quadra_batch_hotspot_benchmark.sh

cat <<'EOF'

Installed HAD Quadra batch propagation hotspot benchmark.

Files added:
  benchmarks/benchmark_had_quadra_batch_hotspots.cpp
  run_had_quadra_batch_hotspot_benchmark.sh

Run:
  ./run_had_quadra_batch_hotspot_benchmark.sh 10

EOF
