#!/usr/bin/env bash
set -euo pipefail

# install_sparse_rw1_hdot_pipeline_timing.sh
#
# Full-pipeline timing benchmark for sparse RW1 Hdot generation.
#
# Buckets:
#   graph_build_ms   : construct AD graph + PropagateAdjoint()
#   seed_ms          : seed K total directions into dotBatch
#   reverse_ms       : PropagateAdjointDirectionalBatch()
#   extract_ms       : read selected Hdot entries from GetAdjointDotBatch()
#   total_ms         : total measured pipeline
#
# Default CXXFLAGS include -g for profiler use.

mkdir -p benchmarks .quadra_patch_backups

if [[ ! -f core/had_quadra.hpp ]]; then
  echo "ERROR: missing core/had_quadra.hpp"
  exit 1
fi

if ! grep -q "PropagateAdjointDirectionalBatch" core/had_quadra.hpp; then
  echo "ERROR: PropagateAdjointDirectionalBatch not found."
  exit 1
fi

cat > benchmarks/benchmark_sparse_rw1_hdot_pipeline_timing.cpp <<'EOF'
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../core/had_quadra.hpp"
#include "../core/laplace/sparse_huu_factorization.hpp"

DECLARE_ADGRAPH()

namespace {

using Clock = std::chrono::steady_clock;
using had::AReal;
using had::Real;

constexpr double kPi = 3.141592653589793238462643383279502884;

double since_ms(const Clock::time_point& a, const Clock::time_point& b) {
    return std::chrono::duration<double, std::milli>(b - a).count();
}

struct SparseRw1Objective {
    int m;

    template <class T>
    T operator()(const std::vector<T>& x) const {
        const T mu = x[0];
        const T log_sigma = x[1];
        const T log_lambda0 = x[2];
        const T log_lambda_rw = x[3];
        const T log_beta = x[4];

        const T inv_sigma2 = exp(T(-2.0) * log_sigma);
        const T lambda0 = exp(log_lambda0);
        const T lambda_rw = exp(log_lambda_rw);
        const T beta = exp(log_beta);

        T nll = T(0.0);

        for (int i = 0; i < m; ++i) {
            const double xd = static_cast<double>(i + 1);
            const T y = T(0.6 + 0.10 * std::sin(0.21 * xd)
                              + 0.07 * std::cos(0.47 * xd));
            const T u = x[5 + i];
            const T resid = y - mu - u;

            nll = nll
                + T(0.5) * resid * resid * inv_sigma2
                + log_sigma
                + T(0.5 * std::log(2.0 * kPi))
                + T(0.5) * lambda0 * u * u
                + beta * exp(u);
        }

        for (int i = 1; i < m; ++i) {
            const T diff = x[5 + i] - x[5 + i - 1];
            nll = nll + T(0.5) * lambda_rw * diff * diff;
        }

        return nll;
    }
};

Eigen::VectorXd make_y(int m) {
    Eigen::VectorXd y(m);
    for (int i = 0; i < m; ++i) {
        const double x = static_cast<double>(i + 1);
        y[i] = 0.6 + 0.10 * std::sin(0.21 * x)
                   + 0.07 * std::cos(0.47 * x);
    }
    return y;
}

Eigen::SparseMatrix<double> Huu_sparse_direct(const Eigen::VectorXd& theta,
                                              const Eigen::VectorXd& u) {
    const int m = static_cast<int>(u.size());
    const double inv_sigma2 = std::exp(-2.0 * theta[1]);
    const double lambda0 = std::exp(theta[2]);
    const double lambda_rw = std::exp(theta[3]);
    const double beta = std::exp(theta[4]);

    std::vector<Eigen::Triplet<double>> triplets;
    triplets.reserve(static_cast<size_t>(3 * m - 2));

    for (int i = 0; i < m; ++i) {
        double diag = inv_sigma2 + lambda0 + beta * std::exp(u[i]);
        if (i > 0) diag += lambda_rw;
        if (i + 1 < m) diag += lambda_rw;
        triplets.emplace_back(i, i, diag);
    }

    for (int i = 1; i < m; ++i) {
        triplets.emplace_back(i, i - 1, -lambda_rw);
        triplets.emplace_back(i - 1, i, -lambda_rw);
    }

    Eigen::SparseMatrix<double> H(m, m);
    H.setFromTriplets(triplets.begin(), triplets.end());
    H.makeCompressed();
    return H;
}

Eigen::VectorXd random_gradient(const Eigen::VectorXd& theta,
                                const Eigen::VectorXd& u) {
    const int m = static_cast<int>(u.size());
    const double mu = theta[0];
    const double inv_sigma2 = std::exp(-2.0 * theta[1]);
    const double lambda0 = std::exp(theta[2]);
    const double lambda_rw = std::exp(theta[3]);
    const double beta = std::exp(theta[4]);
    const Eigen::VectorXd y = make_y(m);

    Eigen::VectorXd g = Eigen::VectorXd::Zero(m);

    for (int i = 0; i < m; ++i) {
        const double resid = y[i] - mu - u[i];
        g[i] += -resid * inv_sigma2
              + lambda0 * u[i]
              + beta * std::exp(u[i]);
    }

    for (int i = 1; i < m; ++i) {
        const double diff = u[i] - u[i - 1];
        g[i] += lambda_rw * diff;
        g[i - 1] -= lambda_rw * diff;
    }

    return g;
}

Eigen::VectorXd solve_uhat(const Eigen::VectorXd& theta, int m) {
    Eigen::VectorXd u = Eigen::VectorXd::Zero(m);

    for (int iter = 0; iter < 80; ++iter) {
        const Eigen::VectorXd g = random_gradient(theta, u);
        Eigen::LDLT<Eigen::MatrixXd> ldlt(Eigen::MatrixXd(Huu_sparse_direct(theta, u)));
        const Eigen::VectorXd step = ldlt.solve(g);
        u -= step;
        if (step.lpNorm<Eigen::Infinity>() < 1.0e-12) break;
    }

    return u;
}

Eigen::VectorXd f_u_theta_column(const Eigen::VectorXd& theta,
                                 const Eigen::VectorXd& uhat,
                                 int theta_index) {
    const int m = static_cast<int>(uhat.size());
    const double inv_sigma2 = std::exp(-2.0 * theta[1]);
    const double lambda0 = std::exp(theta[2]);
    const double lambda_rw = std::exp(theta[3]);
    const double beta = std::exp(theta[4]);
    const Eigen::VectorXd y = make_y(m);

    Eigen::VectorXd col = Eigen::VectorXd::Zero(m);

    if (theta_index == 0) {
        col.array() = inv_sigma2;
        return col;
    }
    if (theta_index == 1) {
        for (int i = 0; i < m; ++i) {
            const double resid = y[i] - theta[0] - uhat[i];
            col[i] = 2.0 * resid * inv_sigma2;
        }
        return col;
    }
    if (theta_index == 2) {
        for (int i = 0; i < m; ++i) col[i] = lambda0 * uhat[i];
        return col;
    }
    if (theta_index == 3) {
        for (int i = 1; i < m; ++i) {
            const double diff = uhat[i] - uhat[i - 1];
            col[i] += lambda_rw * diff;
            col[i - 1] -= lambda_rw * diff;
        }
        return col;
    }
    if (theta_index == 4) {
        for (int i = 0; i < m; ++i) col[i] = beta * std::exp(uhat[i]);
        return col;
    }

    return col;
}

struct Row {
    double graph_build_ms = 0.0;
    double seed_ms = 0.0;
    double reverse_ms = 0.0;
    double extract_ms = 0.0;
    double total_ms = 0.0;
    double checksum = 0.0;
    int vertices = 0;
    int extract_entries = 0;
};

Row run_case(int m, int K) {
    Eigen::VectorXd theta(5);
    theta << 0.55, std::log(0.65), std::log(0.55), std::log(0.90), std::log(0.25);

    const Eigen::VectorXd uhat = solve_uhat(theta, m);
    const Eigen::SparseMatrix<double> H = Huu_sparse_direct(theta, uhat);
    quadra::laplace::SparseHuuFactorization factor(H);

    auto direction_provider = [&](int theta_index) -> Eigen::VectorXd {
        return -factor.solve(f_u_theta_column(theta, uhat, theta_index));
    };

    Row out;

    const auto total_start = Clock::now();

    had::ADGraph graph;
    had::g_ADGraph = &graph;

    std::vector<AReal> x(static_cast<size_t>(5 + m));

    const auto graph_start = Clock::now();

    for (int j = 0; j < 5; ++j) x[static_cast<size_t>(j)] = AReal(theta[j]);
    for (int i = 0; i < m; ++i) x[static_cast<size_t>(5 + i)] = AReal(uhat[i]);

    SparseRw1Objective objective{m};
    AReal f = objective(x);

    had::g_ADGraph->vertices[f.varId].w = 1.0;
    had::PropagateAdjoint();

    const auto graph_end = Clock::now();

    had::ResizeDirectionalBatch(K);

    const auto seed_start = Clock::now();

    for (int k = 0; k < K; ++k) {
        const int theta_index = k % 5;
        const Eigen::VectorXd udir = direction_provider(theta_index);

        for (int j = 0; j < 5; ++j) {
            had::SetARealDotBatch(x[static_cast<size_t>(j)], k,
                                  j == theta_index ? 1.0 : 0.0);
        }
        for (int i = 0; i < m; ++i) {
            had::SetARealDotBatch(x[static_cast<size_t>(5 + i)], k, udir[i]);
        }
    }

    const auto seed_end = Clock::now();

    const auto reverse_start = Clock::now();
    had::PropagateAdjointDirectionalBatch();
    const auto reverse_end = Clock::now();

    const auto extract_start = Clock::now();

    double checksum = 0.0;
    int entries = 0;

    // Extract tridiagonal random-effect Hdot pattern for each direction.
    for (int k = 0; k < K; ++k) {
        for (int i = 0; i < m; ++i) {
            checksum += had::GetAdjointDotBatch(
                x[static_cast<size_t>(5 + i)],
                x[static_cast<size_t>(5 + i)],
                k);
            ++entries;

            if (i > 0) {
                checksum += had::GetAdjointDotBatch(
                    x[static_cast<size_t>(5 + i)],
                    x[static_cast<size_t>(5 + i - 1)],
                    k);
                ++entries;
            }
        }
    }

    const auto extract_end = Clock::now();
    const auto total_end = Clock::now();

    out.graph_build_ms = since_ms(graph_start, graph_end);
    out.seed_ms = since_ms(seed_start, seed_end);
    out.reverse_ms = since_ms(reverse_start, reverse_end);
    out.extract_ms = since_ms(extract_start, extract_end);
    out.total_ms = since_ms(total_start, total_end);
    out.checksum = checksum;
    out.vertices = static_cast<int>(had::g_ADGraph->vertices.size());
    out.extract_entries = entries;

    return out;
}

}  // namespace

int main(int argc, char** argv) {
    int reps = 1;
    if (argc > 1) reps = std::stoi(argv[1]);

    const std::vector<int> m_values = {100, 250, 500};
    const int K = 5;

    std::cout << "Sparse RW1 Hdot pipeline timing benchmark\n";
    std::cout << "reps per case = " << reps << "\n";
    std::cout << "K = " << K << "\n\n";

    std::cout << std::setw(8) << "m"
              << std::setw(12) << "vertices"
              << std::setw(12) << "entries"
              << std::setw(14) << "graph"
              << std::setw(14) << "seed"
              << std::setw(14) << "reverse"
              << std::setw(14) << "extract"
              << std::setw(14) << "total"
              << std::setw(12) << "rev%"
              << std::setw(12) << "ext%"
              << "\n";

    std::cout << std::scientific << std::setprecision(6);

    for (int m : m_values) {
        Row avg;
        for (int r = 0; r < reps; ++r) {
            const Row row = run_case(m, K);
            avg.graph_build_ms += row.graph_build_ms;
            avg.seed_ms += row.seed_ms;
            avg.reverse_ms += row.reverse_ms;
            avg.extract_ms += row.extract_ms;
            avg.total_ms += row.total_ms;
            avg.checksum += row.checksum;
            avg.vertices = row.vertices;
            avg.extract_entries = row.extract_entries;
        }

        avg.graph_build_ms /= reps;
        avg.seed_ms /= reps;
        avg.reverse_ms /= reps;
        avg.extract_ms /= reps;
        avg.total_ms /= reps;

        std::cout << std::setw(8) << m
                  << std::setw(12) << avg.vertices
                  << std::setw(12) << avg.extract_entries
                  << std::setw(14) << avg.graph_build_ms
                  << std::setw(14) << avg.seed_ms
                  << std::setw(14) << avg.reverse_ms
                  << std::setw(14) << avg.extract_ms
                  << std::setw(14) << avg.total_ms
                  << std::setw(12) << (avg.reverse_ms / avg.total_ms)
                  << std::setw(12) << (avg.extract_ms / avg.total_ms)
                  << "\n";
    }

    std::cout << "\nBenchmark complete.\n";
    return 0;
}
EOF

cat > run_sparse_rw1_hdot_pipeline_timing.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O3 -DNDEBUG -g}"
REPS="${1:-1}"

EIGEN_INCLUDE=""
if [[ -d external/Eigen ]]; then
  EIGEN_INCLUDE="-Iexternal/Eigen"
elif [[ -d core/eigen ]]; then
  EIGEN_INCLUDE="-Icore/eigen"
fi

LBFGS_INCLUDE=""
if [[ -d external/LBFGSpp/include ]]; then
  LBFGS_INCLUDE="-Iexternal/LBFGSpp/include"
elif [[ -d external/LBFGSpp ]]; then
  LBFGS_INCLUDE="-Iexternal/LBFGSpp"
fi

mkdir -p build/benchmarks

set -x
"${CXX}" ${CXXFLAGS} ${EIGEN_INCLUDE} ${LBFGS_INCLUDE} -I. \
  benchmarks/benchmark_sparse_rw1_hdot_pipeline_timing.cpp \
  -o build/benchmarks/benchmark_sparse_rw1_hdot_pipeline_timing

./build/benchmarks/benchmark_sparse_rw1_hdot_pipeline_timing "${REPS}"
EOF

chmod +x run_sparse_rw1_hdot_pipeline_timing.sh

cat <<'EOF'

Installed sparse RW1 Hdot pipeline timing benchmark.

Files added:
  benchmarks/benchmark_sparse_rw1_hdot_pipeline_timing.cpp
  run_sparse_rw1_hdot_pipeline_timing.sh

Run:
  ./run_had_quadra_directional_batch_propagation_test.sh
  ./run_sparse_rw1_hdot_pipeline_timing.sh 1

Profiler-friendly binary:
  build/benchmarks/benchmark_sparse_rw1_hdot_pipeline_timing

EOF
