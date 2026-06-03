#!/usr/bin/env bash
set -euo pipefail

# install_fixed_theta_quadra_tmb_runtime_benchmark_v1.sh
#
# Adds fixed-theta runtime benchmarks:
#
#   1. Quadra dense finite-difference Laplace repeated evals
#   2. TMB fixed-theta repeated obj$fn() evals
#
# This benchmark intentionally does NOT optimize fixed effects. It compares
# repeated marginal objective evaluation for the same model/parameters.

mkdir -p examples/state_space_surplus_production examples/tmb_state_space_surplus

if [[ ! -f examples/state_space_surplus_production/state_space_surplus_production.hpp ]]; then
  echo "ERROR: missing state-space surplus header"
  exit 1
fi

cat > examples/state_space_surplus_production/benchmark_dense_laplace_fixed_theta.cpp <<'EOF'
#include "state_space_surplus_production.hpp"

#include <Eigen/Dense>
#include <LBFGS.h>

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <vector>

namespace ss = quadra_examples::state_space_surplus_production;
using Clock = std::chrono::high_resolution_clock;

namespace {

std::vector<double> to_std_vector(const Eigen::VectorXd& x) {
  std::vector<double> out(static_cast<std::size_t>(x.size()));
  for (int i = 0; i < x.size(); ++i) out[static_cast<std::size_t>(i)] = x[i];
  return out;
}

Eigen::VectorXd to_eigen_vector(const std::vector<double>& x) {
  Eigen::VectorXd out(static_cast<int>(x.size()));
  for (std::size_t i = 0; i < x.size(); ++i) out[static_cast<int>(i)] = x[i];
  return out;
}

double ms_between(const Clock::time_point& a, const Clock::time_point& b) {
  return std::chrono::duration<double, std::milli>(b - a).count();
}

double joint_u(const ss::Data& data, const ss::Parameters& par, const Eigen::VectorXd& u) {
  return ss::joint_objective(data, par, to_std_vector(u));
}

Eigen::VectorXd fd_grad_u(const ss::Data& data, const ss::Parameters& par, const Eigen::VectorXd& u) {
  Eigen::VectorXd grad(u.size());
  for (int i = 0; i < u.size(); ++i) {
    const double h = 1e-5 * (1.0 + std::abs(u[i]));
    Eigen::VectorXd up = u, um = u;
    up[i] += h;
    um[i] -= h;
    grad[i] = (joint_u(data, par, up) - joint_u(data, par, um)) / (2.0 * h);
  }
  return grad;
}

class ObjU {
 public:
  ObjU(const ss::Data& data, const ss::Parameters& par) : data_(data), par_(par) {}
  double operator()(const Eigen::VectorXd& u, Eigen::VectorXd& grad) {
    const double f = joint_u(data_, par_, u);
    grad = fd_grad_u(data_, par_, u);
    return f;
  }
 private:
  const ss::Data& data_;
  const ss::Parameters& par_;
};

Eigen::VectorXd optimize_u(const ss::Data& data, const ss::Parameters& par) {
  Eigen::VectorXd u = to_eigen_vector(ss::zero_random_effects(data));

  LBFGSpp::LBFGSParam<double> param;
  param.epsilon = 1e-7;
  param.max_iterations = 500;
  param.max_linesearch = 100;
  param.m = 8;
  param.ftol = 1e-4;
  param.wolfe = 0.9;
  param.min_step = 1e-20;
  param.max_step = 1.0;

  LBFGSpp::LBFGSSolver<double> solver(param);
  ObjU obj(data, par);
  double f = 0.0;

  try {
    solver.minimize(obj, u, f);
  } catch (...) {
    // Accept the same benign finite-difference line-search termination used
    // by the example if the gradient is already small.
    const double gnorm = fd_grad_u(data, par, u).norm();
    if (!(std::isfinite(joint_u(data, par, u)) && gnorm < 1e-3)) {
      throw;
    }
  }

  return u;
}

Eigen::MatrixXd fd_hessian_uu(const ss::Data& data, const ss::Parameters& par, const Eigen::VectorXd& u) {
  const int n = static_cast<int>(u.size());
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(n, n);
  const double f0 = joint_u(data, par, u);

  for (int i = 0; i < n; ++i) {
    const double hi = 1e-4 * (1.0 + std::abs(u[i]));
    Eigen::VectorXd up = u, um = u;
    up[i] += hi;
    um[i] -= hi;
    H(i, i) = (joint_u(data, par, up) - 2.0 * f0 + joint_u(data, par, um)) / (hi * hi);

    for (int j = i + 1; j < n; ++j) {
      const double hj = 1e-4 * (1.0 + std::abs(u[j]));
      Eigen::VectorXd upp = u, upm = u, ump = u, umm = u;
      upp[i] += hi; upp[j] += hj;
      upm[i] += hi; upm[j] -= hj;
      ump[i] -= hi; ump[j] += hj;
      umm[i] -= hi; umm[j] -= hj;
      const double hij =
          (joint_u(data, par, upp) - joint_u(data, par, upm) -
           joint_u(data, par, ump) + joint_u(data, par, umm)) / (4.0 * hi * hj);
      H(i, j) = hij;
      H(j, i) = hij;
    }
  }

  return H;
}

double dense_laplace_eval(const ss::Data& data, const ss::Parameters& par) {
  const Eigen::VectorXd uhat = optimize_u(data, par);
  const double joint = joint_u(data, par, uhat);
  const Eigen::MatrixXd H = fd_hessian_uu(data, par, uhat);

  Eigen::LLT<Eigen::MatrixXd> llt(H);
  if (llt.info() != Eigen::Success) {
    throw std::runtime_error("Huu not SPD");
  }

  const auto& L = llt.matrixL();
  double logdet = 0.0;
  for (int i = 0; i < H.rows(); ++i) {
    logdet += 2.0 * std::log(L(i, i));
  }

  const double n_u = static_cast<double>(uhat.size());
  return joint + 0.5 * logdet - 0.5 * n_u * std::log(2.0 * M_PI);
}

}  // namespace

int main(int argc, char** argv) {
  int reps = 20;
  if (argc > 1) reps = std::stoi(argv[1]);

  const ss::Data data = ss::make_demo_data();
  const ss::Parameters par = ss::make_demo_parameters();

  double last = 0.0;

  // warmup
  last = dense_laplace_eval(data, par);

  const auto t0 = Clock::now();
  for (int i = 0; i < reps; ++i) {
    last = dense_laplace_eval(data, par);
  }
  const auto t1 = Clock::now();

  const double total_ms = ms_between(t0, t1);
  const double avg_ms = total_ms / static_cast<double>(reps);

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "Quadra dense finite-difference Laplace fixed-theta benchmark\n";
  std::cout << "reps = " << reps << "\n";
  std::cout << "objective = " << last << "\n";
  std::cout << "total_ms = " << total_ms << "\n";
  std::cout << "avg_ms = " << avg_ms << "\n";

  return 0;
}
EOF

cat > run_quadra_dense_laplace_fixed_theta_benchmark.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-20}"
CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/examples

set -x
"${CXX}" ${CXXFLAGS} \
  -Iexternal/Eigen \
  -Iexternal/LBFGSpp/include \
  -Iexamples/state_space_surplus_production \
  -Iexamples/surplus_production \
  examples/state_space_surplus_production/benchmark_dense_laplace_fixed_theta.cpp \
  -o build/examples/benchmark_dense_laplace_fixed_theta

./build/examples/benchmark_dense_laplace_fixed_theta "$REPS"
EOF

chmod +x run_quadra_dense_laplace_fixed_theta_benchmark.sh

cat > examples/tmb_state_space_surplus/benchmark_fixed_theta_tmb.R <<'EOF'
#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[[1]]) else 20L

cat("TMB fixed-theta state-space surplus benchmark\n")
cat("============================================\n")

if (!requireNamespace("TMB", quietly = TRUE)) {
  cat("TMB is not installed. Skipping.\n")
  quit(status = 0)
}

library(TMB)

template <- file.path("examples", "tmb_state_space_surplus", "state_space_surplus_tmb.cpp")
dynlib_name <- "state_space_surplus_tmb"

catch_observed <- c(
  80, 88, 95, 105, 115, 125, 130, 128,
  120, 110, 100, 90, 82, 78, 75
)

index_observed <- c(
  1.55, 1.50, 1.43, 1.34, 1.22, 1.10, 0.98, 0.88,
  0.81, 0.78, 0.77, 0.80, 0.84, 0.88, 0.92
)

parameters <- list(
  log_r = log(0.5),
  log_K = log(700.0),
  log_q = log(0.0024),
  log_sigma_process = log(0.15),
  log_sigma_index = log(0.10),
  logit_B0_frac = log(0.90 / 0.10),
  u = rep(0, length(catch_observed) - 1)
)

data <- list(catch_observed = catch_observed, index_observed = index_observed)

if (!file.exists(TMB::dynlib(file.path("examples", "tmb_state_space_surplus", dynlib_name)))) {
  cat("Compiling TMB template...\n")
  TMB::compile(template, flags = "-O2")
}

dyn.load(TMB::dynlib(file.path("examples", "tmb_state_space_surplus", dynlib_name)))

obj <- TMB::MakeADFun(
  data = data,
  parameters = parameters,
  random = "u",
  DLL = dynlib_name,
  silent = TRUE
)

# Warmup.
last <- obj$fn()

gc()
t0 <- proc.time()
for (i in seq_len(reps)) {
  last <- obj$fn()
}
t1 <- proc.time()

elapsed <- as.numeric((t1 - t0)[["elapsed"]])
avg_ms <- elapsed * 1000 / reps

cat("reps =", reps, "\n")
cat("objective =", format(last, digits = 12), "\n")
cat("total_ms =", format(elapsed * 1000, digits = 12), "\n")
cat("avg_ms =", format(avg_ms, digits = 12), "\n")
cat("expected Quadra dense FD objective ~= -10.642184\n")
EOF

chmod +x examples/tmb_state_space_surplus/benchmark_fixed_theta_tmb.R

cat > run_tmb_fixed_theta_benchmark.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-20}"
Rscript examples/tmb_state_space_surplus/benchmark_fixed_theta_tmb.R "$REPS"
EOF

chmod +x run_tmb_fixed_theta_benchmark.sh

cat > run_quadra_vs_tmb_fixed_theta_benchmark.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-20}"

echo "== Quadra dense finite-difference Laplace =="
./run_quadra_dense_laplace_fixed_theta_benchmark.sh "$REPS"

echo
echo "== TMB AD/Laplace =="
./run_tmb_fixed_theta_benchmark.sh "$REPS"
EOF

chmod +x run_quadra_vs_tmb_fixed_theta_benchmark.sh

cat <<'EOF'

Installed fixed-theta Quadra vs TMB runtime benchmark.

Added:
  examples/state_space_surplus_production/benchmark_dense_laplace_fixed_theta.cpp
  run_quadra_dense_laplace_fixed_theta_benchmark.sh
  examples/tmb_state_space_surplus/benchmark_fixed_theta_tmb.R
  run_tmb_fixed_theta_benchmark.sh
  run_quadra_vs_tmb_fixed_theta_benchmark.sh

Run:
  ./run_quadra_vs_tmb_fixed_theta_benchmark.sh 20

Note:
  Quadra path is still dense finite-difference Laplace.
  TMB path is AD/Laplace.
  This is a baseline comparison, not yet a fair final Quadra exact-derivative benchmark.

EOF
