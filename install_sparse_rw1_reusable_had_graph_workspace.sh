#!/usr/bin/env bash
set -euo pipefail

# install_sparse_rw1_reusable_had_graph_workspace.sh
#
# Adds a benchmark scaffold for graph reuse:
#
#   build AD graph once
#   reuse graph object and AReal handles across repeated Hdot evaluations
#
# This is intentionally a benchmark/scaffold first, not a production API.
# It compares:
#
#   rebuild-each-time pipeline
#   vs
#   reuse-workspace pipeline
#
# for sparse RW1 Hdot timing.

mkdir -p benchmarks .quadra_patch_backups

if [[ ! -f core/had_quadra.hpp ]]; then
  echo "ERROR: missing core/had_quadra.hpp"
  exit 1
fi

if ! grep -q "PropagateAdjointDirectionalBatch" core/had_quadra.hpp; then
  echo "ERROR: PropagateAdjointDirectionalBatch not found."
  exit 1
fi

cat > benchmarks/benchmark_sparse_rw1_reusable_had_graph_workspace.cpp <<'EOF'
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

double ms_between(const Clock::time_point& a, const Clock::time_point& b) {
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

struct ReusableRw1HadWorkspace {
    int m = 0;
    int K = 0;
    had::ADGraph graph;
    std::vector<AReal> x;
    AReal f;

    ReusableRw1HadWorkspace(int m_, int K_,
                            const Eigen::VectorXd& theta,
                            const Eigen::VectorXd& uhat)
        : m(m_), K(K_), x(static_cast<size_t>(5 + m_)) {
        build(theta, uhat);
    }

    void build(const Eigen::VectorXd& theta,
               const Eigen::VectorXd& uhat) {
        had::g_ADGraph = &graph;

        for (int j = 0; j < 5; ++j) {
            x[static_cast<size_t>(j)] = AReal(theta[j]);
        }
        for (int i = 0; i < m; ++i) {
            x[static_cast<size_t>(5 + i)] = AReal(uhat[i]);
        }

        SparseRw1Objective objective{m};
        f = objective(x);

        had::ResizeDirectionalBatch(K);
    }

    void reset_adjoint_and_second_order_storage() {
        had::g_ADGraph = &graph;

        for (auto& v : graph.vertices) {
            v.w = Real(0.0);
            v.wDot = Real(0.0);
            v.soW = v.soW;      // keep structural/local weights
            v.soWDot = Real(0.0);
        }

        if (graph.soEdges.size() < graph.vertices.size()) {
            graph.soEdges.resize(graph.vertices.size());
        } else {
            for (auto& tree : graph.soEdges) tree.Clear();
        }

        if (graph.soEdgesDot.size() < graph.vertices.size()) {
            graph.soEdgesDot.resize(graph.vertices.size());
        } else {
            for (auto& tree : graph.soEdgesDot) tree.Clear();
        }

        graph.selfSoEdges.assign(graph.vertices.size(), Real(0.0));
        graph.selfSoEdgesDot.assign(graph.vertices.size(), Real(0.0));
    }

    void propagate_base_adjoint() {
        had::g_ADGraph = &graph;
        reset_adjoint_and_second_order_storage();
        graph.vertices[f.varId].w = Real(1.0);
        had::PropagateAdjoint();
    }

    void seed_directions(const Eigen::VectorXd& theta,
                         const Eigen::VectorXd& uhat) {
        (void)theta;

        const Eigen::SparseMatrix<double> H = Huu_sparse_direct(theta, uhat);
        quadra::laplace::SparseHuuFactorization factor(H);

        auto direction_provider = [&](int theta_index) -> Eigen::VectorXd {
            return -factor.solve(f_u_theta_column(theta, uhat, theta_index));
        };

        had::g_ADGraph = &graph;
        had::ResizeDirectionalBatch(K);

        for (int k = 0; k < K; ++k) {
            const int theta_index = k % 5;
            const Eigen::VectorXd udir = direction_provider(theta_index);

            for (int j = 0; j < 5; ++j) {
                had::SetARealDotBatch(x[static_cast<size_t>(j)], k,
                                      j == theta_index ? 1.0 : 0.0);
            }
            for (int i = 0; i < m; ++i) {
                had::SetARealDotBatch(x[static_cast<size_t>(5 + i)], k,
                                      udir[i]);
            }
        }
    }

    double extract_checksum() {
        had::g_ADGraph = &graph;

        double checksum = 0.0;

        for (int k = 0; k < K; ++k) {
            for (int i = 0; i < m; ++i) {
                checksum += had::GetAdjointDotBatch(
                    x[static_cast<size_t>(5 + i)],
                    x[static_cast<size_t>(5 + i)],
                    k);

                if (i > 0) {
                    checksum += had::GetAdjointDotBatch(
                        x[static_cast<size_t>(5 + i)],
                        x[static_cast<size_t>(5 + i - 1)],
                        k);
                }
            }
        }

        return checksum;
    }
};

struct Row {
    double rebuild_ms = 0.0;
    double reuse_ms = 0.0;
    double speedup = 0.0;
    double checksum = 0.0;
    int vertices = 0;
};

double run_rebuild_once(int m, int K,
                        const Eigen::VectorXd& theta,
                        const Eigen::VectorXd& uhat) {
    had::ADGraph graph;
    had::g_ADGraph = &graph;

    std::vector<AReal> x(static_cast<size_t>(5 + m));

    for (int j = 0; j < 5; ++j) x[static_cast<size_t>(j)] = AReal(theta[j]);
    for (int i = 0; i < m; ++i) x[static_cast<size_t>(5 + i)] = AReal(uhat[i]);

    SparseRw1Objective objective{m};
    AReal f = objective(x);

    had::g_ADGraph->vertices[f.varId].w = Real(1.0);
    had::PropagateAdjoint();

    const Eigen::SparseMatrix<double> H = Huu_sparse_direct(theta, uhat);
    quadra::laplace::SparseHuuFactorization factor(H);

    auto direction_provider = [&](int theta_index) -> Eigen::VectorXd {
        return -factor.solve(f_u_theta_column(theta, uhat, theta_index));
    };

    had::ResizeDirectionalBatch(K);

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

    had::PropagateAdjointDirectionalBatch();

    double checksum = 0.0;
    for (int k = 0; k < K; ++k) {
        for (int i = 0; i < m; ++i) {
            checksum += had::GetAdjointDotBatch(
                x[static_cast<size_t>(5 + i)],
                x[static_cast<size_t>(5 + i)],
                k);
            if (i > 0) {
                checksum += had::GetAdjointDotBatch(
                    x[static_cast<size_t>(5 + i)],
                    x[static_cast<size_t>(5 + i - 1)],
                    k);
            }
        }
    }

    return checksum;
}

Row run_case(int m, int K, int reps) {
    Eigen::VectorXd theta(5);
    theta << 0.55, std::log(0.65), std::log(0.55), std::log(0.90), std::log(0.25);
    const Eigen::VectorXd uhat = solve_uhat(theta, m);

    Row out;

    auto t0 = Clock::now();
    double rebuild_checksum = 0.0;
    for (int r = 0; r < reps; ++r) {
        rebuild_checksum += run_rebuild_once(m, K, theta, uhat);
    }
    auto t1 = Clock::now();

    ReusableRw1HadWorkspace workspace(m, K, theta, uhat);
    workspace.propagate_base_adjoint();

    auto t2 = Clock::now();
    double reuse_checksum = 0.0;
    for (int r = 0; r < reps; ++r) {
        workspace.propagate_base_adjoint();
        workspace.seed_directions(theta, uhat);
        had::PropagateAdjointDirectionalBatch();
        reuse_checksum += workspace.extract_checksum();
    }
    auto t3 = Clock::now();

    out.rebuild_ms = ms_between(t0, t1) / static_cast<double>(reps);
    out.reuse_ms = ms_between(t2, t3) / static_cast<double>(reps);
    out.speedup = out.rebuild_ms / out.reuse_ms;
    out.checksum = std::abs(rebuild_checksum - reuse_checksum);
    out.vertices = static_cast<int>(workspace.graph.vertices.size());

    return out;
}

}  // namespace

int main(int argc, char** argv) {
    int reps = 10;
    if (argc > 1) reps = std::stoi(argv[1]);

    const std::vector<int> m_values = {100, 250, 500};
    const int K = 5;

    std::cout << "Sparse RW1 reusable HAD graph workspace benchmark\n";
    std::cout << "reps per case = " << reps << "\n";
    std::cout << "K = " << K << "\n\n";

    std::cout << std::setw(8) << "m"
              << std::setw(12) << "vertices"
              << std::setw(16) << "rebuild ms"
              << std::setw(16) << "reuse ms"
              << std::setw(14) << "speedup"
              << std::setw(16) << "checksum diff"
              << "\n";

    std::cout << std::scientific << std::setprecision(6);

    for (int m : m_values) {
        const Row r = run_case(m, K, reps);

        std::cout << std::setw(8) << m
                  << std::setw(12) << r.vertices
                  << std::setw(16) << r.rebuild_ms
                  << std::setw(16) << r.reuse_ms
                  << std::setw(14) << r.speedup
                  << std::setw(16) << r.checksum
                  << "\n";
    }

    std::cout << "\nBenchmark complete.\n";
    return 0;
}
EOF

cat > run_sparse_rw1_reusable_had_graph_workspace.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O3 -DNDEBUG -g}"
REPS="${1:-10}"

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
  benchmarks/benchmark_sparse_rw1_reusable_had_graph_workspace.cpp \
  -o build/benchmarks/benchmark_sparse_rw1_reusable_had_graph_workspace

./build/benchmarks/benchmark_sparse_rw1_reusable_had_graph_workspace "${REPS}"
EOF

chmod +x run_sparse_rw1_reusable_had_graph_workspace.sh

cat <<'EOF'

Installed sparse RW1 reusable HAD graph workspace benchmark.

Files added:
  benchmarks/benchmark_sparse_rw1_reusable_had_graph_workspace.cpp
  run_sparse_rw1_reusable_had_graph_workspace.sh

Run:
  ./run_had_quadra_directional_batch_propagation_test.sh
  ./run_sparse_rw1_reusable_had_graph_workspace.sh 10

EOF
