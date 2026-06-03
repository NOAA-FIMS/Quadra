#!/usr/bin/env bash
set -euo pipefail

# install_sparse_rw1_flat_accumulator_trace_check_v1.sh
#
# Adds a focused integration check for the flat accumulator:
#
#   BTree-backed exact-gradient reuse trace terms
#        vs
#   BatchDirectionalFlatAccumulator RW1 trace terms
#
# This is the bridge between the microbenchmark and a future specialized
# PropagateAdjointDirectionalBatchFlat() backend.
#
# It does not alter production reverse propagation.

mkdir -p benchmarks .quadra_patch_backups

target="benchmarks/benchmark_sparse_rw1_flat_accumulator_trace_check.cpp"
if [[ -f "$target" ]]; then
  cp "$target" ".quadra_patch_backups/benchmark_sparse_rw1_flat_accumulator_trace_check.cpp.$(date +%Y%m%d_%H%M%S).bak"
fi

cat > "$target" <<'EOF'
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "../core/had/batch_directional_flat_accumulator.hpp"
#include "../core/had_quadra.hpp"
#include "../core/laplace/exact_gradient_workspace.hpp"
#include "../core/laplace/sparse_huu_factorization.hpp"
#include "../core/laplace/structure_aware_rw1_hdot.hpp"

DECLARE_ADGRAPH()

namespace {

using Clock = std::chrono::steady_clock;
using had::AReal;

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
      const T y =
          T(0.6 + 0.10 * std::sin(0.21 * xd) + 0.07 * std::cos(0.47 * xd));
      const T u = x[5 + i];
      const T resid = y - mu - u;

      nll = nll + T(0.5) * resid * resid * inv_sigma2 + log_sigma +
            T(0.5 * std::log(2.0 * 3.14159265358979323846)) +
            T(0.5) * lambda0 * u * u + beta * exp(u);
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
    y[i] = 0.6 + 0.10 * std::sin(0.21 * x) + 0.07 * std::cos(0.47 * x);
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
  triplets.reserve(static_cast<std::size_t>(3 * m - 2));

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
    g[i] += -resid * inv_sigma2 + lambda0 * u[i] + beta * std::exp(u[i]);
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

    if (step.lpNorm<Eigen::Infinity>() < 1.0e-12) {
      break;
    }
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

struct SelectedTridiagonalInverse {
  Eigen::VectorXd diag;
  Eigen::VectorXd subdiag;

  double operator()(int row, int col) const {
    if (row == col) return diag[row];

    const int hi = std::max(row, col);
    const int lo = std::min(row, col);

    if (hi == lo + 1) return subdiag[hi - 1];

    return 0.0;
  }
};

SelectedTridiagonalInverse compute_selected_inverse(
    Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>>& ldlt, int m) {
  SelectedTridiagonalInverse out;
  out.diag = Eigen::VectorXd::Zero(m);
  out.subdiag = Eigen::VectorXd::Zero(std::max(0, m - 1));

  Eigen::VectorXd rhs = Eigen::VectorXd::Zero(m);

  for (int j = 0; j < m; ++j) {
    rhs.setZero();
    rhs[j] = 1.0;
    const Eigen::VectorXd col = ldlt.solve(rhs);
    out.diag[j] = col[j];
    if (j > 0) out.subdiag[j - 1] = col[j - 1];
  }

  return out;
}

struct BTreeTraceResult {
  Eigen::VectorXd traces;
  double ms = 0.0;
  std::uint64_t queries = 0;
  std::uint64_t pushdots = 0;
  std::uint64_t inserts = 0;
};

BTreeTraceResult btree_hdot_traces(int m,
                                   int K,
                                   const Eigen::VectorXd& theta,
                                   const Eigen::VectorXd& uhat,
                                   const SelectedTridiagonalInverse& selected,
                                   quadra::laplace::SparseHuuFactorization& factor) {
  quadra::laplace::ExactGradientWorkspace workspace;
  std::vector<AReal> theta_vars(5);
  std::vector<AReal> random_vars(static_cast<std::size_t>(m));

  workspace.Build(
      [&]() {
        for (int j = 0; j < 5; ++j) theta_vars[j] = AReal(theta[j]);
        for (int i = 0; i < m; ++i) random_vars[static_cast<std::size_t>(i)] = AReal(uhat[i]);

        std::vector<AReal> x;
        x.reserve(static_cast<std::size_t>(5 + m));
        for (auto& v : theta_vars) x.push_back(v);
        for (auto& v : random_vars) x.push_back(v);

        return SparseRw1Objective{m}(x);
      },
      &theta_vars,
      &random_vars);

  workspace.ResizeDirectionalBatch(static_cast<std::size_t>(K));

  const auto t0 = Clock::now();

  workspace.PropagateBaseAdjoint();

  workspace.SeedTotalDirections(
      static_cast<std::size_t>(K),
      [&](std::size_t k,
          Eigen::VectorXd& theta_direction,
          Eigen::VectorXd& random_direction) {
        const int theta_index = static_cast<int>(k % 5);
        theta_direction = Eigen::VectorXd::Zero(5);
        theta_direction[theta_index] = 1.0;

        random_direction =
            -factor.solve(f_u_theta_column(theta, uhat, theta_index));
      });

  had::ResetBatchDirectionalCounters();
  workspace.PropagateDirectionalBatch();

  Eigen::VectorXd traces = Eigen::VectorXd::Zero(K);

  for (int k = 0; k < K; ++k) {
    double trace = 0.0;

    for (int i = 0; i < m; ++i) {
      trace += selected(i, i) *
               had::GetAdjointDotBatch(random_vars[static_cast<std::size_t>(i)],
                                        random_vars[static_cast<std::size_t>(i)],
                                        k);

      if (i > 0) {
        trace += 2.0 * selected(i, i - 1) *
                 had::GetAdjointDotBatch(random_vars[static_cast<std::size_t>(i)],
                                          random_vars[static_cast<std::size_t>(i - 1)],
                                          k);
      }
    }

    traces[k] = trace;
  }

  const auto t1 = Clock::now();

  BTreeTraceResult out;
  out.traces = traces;
  out.ms = ms_between(t0, t1);
  out.queries = had::g_batch_query_count;
  out.pushdots = had::g_batch_pushdot_count;
  out.inserts = had::g_batch_insert_count;
  return out;
}

struct FlatTraceResult {
  Eigen::VectorXd traces;
  double ms = 0.0;
};

FlatTraceResult flat_accumulator_hdot_traces(
    int m,
    int K,
    const Eigen::VectorXd& theta,
    const Eigen::VectorXd& uhat,
    const SelectedTridiagonalInverse& selected,
    quadra::laplace::SparseHuuFactorization& factor) {
  const auto t0 = Clock::now();

  // 2*m-1 slots: diag entries and subdiag entries.
  had::BatchDirectionalFlatAccumulator accumulator(
      static_cast<std::size_t>(K),
      static_cast<std::size_t>(2 * m - 1));

  accumulator.Clear();

  const double lambda0 = std::exp(theta[2]);
  const double beta = std::exp(theta[4]);

  for (int k = 0; k < K; ++k) {
    const int theta_index = k % 5;

    Eigen::VectorXd theta_direction = Eigen::VectorXd::Zero(5);
    theta_direction[theta_index] = 1.0;

    const Eigen::VectorXd u_direction =
        -factor.solve(f_u_theta_column(theta, uhat, theta_index));

    for (int i = 0; i < m; ++i) {
      double hdot_diag = 0.0;

      if (theta_direction[2] != 0.0) {
        hdot_diag += lambda0 * theta_direction[2];
      }

      hdot_diag += beta * std::exp(uhat[i]) *
                   (theta_direction[4] + u_direction[i]);

      accumulator.Add(static_cast<std::size_t>(k),
                      static_cast<std::size_t>(i),
                      hdot_diag);

      if (i > 0) {
        const int subdiag_slot = m + i - 1;
        const double hdot_subdiag = 0.0;
        accumulator.Add(static_cast<std::size_t>(k),
                        static_cast<std::size_t>(subdiag_slot),
                        hdot_subdiag);
      }
    }
  }

  Eigen::VectorXd traces = Eigen::VectorXd::Zero(K);

  for (int k = 0; k < K; ++k) {
    double trace = 0.0;

    for (int i = 0; i < m; ++i) {
      trace += selected(i, i) *
               accumulator(static_cast<std::size_t>(k),
                           static_cast<std::size_t>(i));

      if (i > 0) {
        const int subdiag_slot = m + i - 1;
        trace += 2.0 * selected(i, i - 1) *
                 accumulator(static_cast<std::size_t>(k),
                             static_cast<std::size_t>(subdiag_slot));
      }
    }

    traces[k] = trace;
  }

  const auto t1 = Clock::now();

  FlatTraceResult out;
  out.traces = traces;
  out.ms = ms_between(t0, t1);
  return out;
}

struct Row {
  int m = 0;
  double btree_ms = 0.0;
  double flat_ms = 0.0;
  double speedup = 0.0;
  double trace_diff = 0.0;
  std::uint64_t queries = 0;
  std::uint64_t pushdots = 0;
  std::uint64_t inserts = 0;
};

Row run_case(int m, int K, int reps) {
  Eigen::VectorXd theta(5);
  theta << 0.55, std::log(0.65), std::log(0.55), std::log(0.90), std::log(0.25);

  const Eigen::VectorXd uhat = solve_uhat(theta, m);
  const Eigen::SparseMatrix<double> H = Huu_sparse_direct(theta, uhat);

  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt(H);
  auto selected = compute_selected_inverse(ldlt, m);

  quadra::laplace::SparseHuuFactorization factor(H);

  Row row;
  row.m = m;

  for (int r = 0; r < reps; ++r) {
    const auto btree = btree_hdot_traces(m, K, theta, uhat, selected, factor);
    const auto flat = flat_accumulator_hdot_traces(m, K, theta, uhat, selected, factor);

    row.btree_ms += btree.ms;
    row.flat_ms += flat.ms;
    row.trace_diff = std::max(row.trace_diff,
                              (btree.traces - flat.traces).cwiseAbs().maxCoeff());
    row.queries = btree.queries;
    row.pushdots = btree.pushdots;
    row.inserts = btree.inserts;
  }

  row.btree_ms /= static_cast<double>(reps);
  row.flat_ms /= static_cast<double>(reps);
  row.speedup = row.btree_ms / row.flat_ms;

  return row;
}

}  // namespace

int main(int argc, char** argv) {
  int reps = 10;
  if (argc > 1) reps = std::stoi(argv[1]);

  const int K = 5;

  std::cout << "Sparse RW1 flat accumulator trace check\n";
  std::cout << "reps per case = " << reps << "\n";
  std::cout << "K = " << K << "\n\n";

  std::cout << std::setw(8) << "m"
            << std::setw(14) << "btree ms"
            << std::setw(14) << "flat ms"
            << std::setw(14) << "speedup"
            << std::setw(14) << "trace diff"
            << std::setw(14) << "queries"
            << std::setw(14) << "pushdots"
            << std::setw(14) << "inserts"
            << "\n";

  std::cout << std::scientific << std::setprecision(6);

  for (int m : {100, 250, 500}) {
    const Row row = run_case(m, K, reps);

    std::cout << std::setw(8) << row.m
              << std::setw(14) << row.btree_ms
              << std::setw(14) << row.flat_ms
              << std::setw(14) << row.speedup
              << std::setw(14) << row.trace_diff
              << std::setw(14) << row.queries
              << std::setw(14) << row.pushdots
              << std::setw(14) << row.inserts
              << "\n";
  }

  std::cout << "\nBenchmark complete.\n";
  return 0;
}
EOF

cat > run_sparse_rw1_flat_accumulator_trace_check.sh <<'EOF'
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
  benchmarks/benchmark_sparse_rw1_flat_accumulator_trace_check.cpp \
  -o build/benchmarks/benchmark_sparse_rw1_flat_accumulator_trace_check

./build/benchmarks/benchmark_sparse_rw1_flat_accumulator_trace_check "${REPS}"
EOF

chmod +x run_sparse_rw1_flat_accumulator_trace_check.sh

cat <<'EOF'

Installed sparse RW1 flat accumulator trace check.

Run:
  ./run_sparse_rw1_flat_accumulator_trace_check.sh 10

Expected:
  trace diff = 0 or near machine precision
  flat ms much smaller than btree ms

EOF
