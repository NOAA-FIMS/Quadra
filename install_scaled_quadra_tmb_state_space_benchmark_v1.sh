#!/usr/bin/env bash
set -euo pipefail

mkdir -p examples/state_space_surplus_production examples/tmb_state_space_surplus

cat > examples/state_space_surplus_production/benchmark_latent_tridiagonal_scaled.cpp <<'EOF'
#include "state_space_surplus_production.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <LBFGS.h>

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

namespace ss = quadra_examples::state_space_surplus_production;
using Clock = std::chrono::high_resolution_clock;

namespace {

double ms_between(const Clock::time_point& a, const Clock::time_point& b) {
  return std::chrono::duration<double, std::milli>(b - a).count();
}

std::vector<int> parse_lengths(const std::string& s) {
  std::vector<int> out;
  std::stringstream ss_in(s);
  std::string item;
  while (std::getline(ss_in, item, ',')) {
    if (!item.empty()) out.push_back(std::stoi(item));
  }
  return out;
}

double inv_logit(const double x) {
  return 1.0 / (1.0 + std::exp(-x));
}

double normal_nll(const double residual, const double sigma) {
  const double z = residual / sigma;
  return std::log(sigma) + 0.5 * std::log(2.0 * M_PI) + 0.5 * z * z;
}

ss::Data make_scaled_data(const int n) {
  ss::Data data;
  data.catch_observed.resize(static_cast<std::size_t>(n));
  data.index_observed.resize(static_cast<std::size_t>(n));

  const double r = 0.5;
  const double K = 700.0;
  const double q = 0.0024;
  double B = 0.90 * K;

  for (int t = 0; t < n; ++t) {
    const double seasonal = std::sin(2.0 * M_PI * static_cast<double>(t) / 17.0);
    const double trend = 1.0 + 0.10 * std::sin(2.0 * M_PI * static_cast<double>(t) / 53.0);
    const double C = 88.0 * trend + 18.0 * seasonal;

    data.catch_observed[static_cast<std::size_t>(t)] = std::max(40.0, C);

    const double obs_error =
        0.05 * std::sin(2.0 * M_PI * static_cast<double>(t) / 11.0) +
        0.025 * std::cos(2.0 * M_PI * static_cast<double>(t) / 7.0);

    data.index_observed[static_cast<std::size_t>(t)] = q * B * std::exp(obs_error);

    if (t < n - 1) {
      const double production = r * B * (1.0 - B / K);
      B = std::max(B + production - data.catch_observed[static_cast<std::size_t>(t)], 1e-9);
    }
  }

  return data;
}

struct Transformed {
  double r;
  double K;
  double q;
  double sigma_process;
  double sigma_index;
  double B0_frac;
  double B0;
};

Transformed transform(const ss::Parameters& par) {
  Transformed out;
  out.r = std::exp(par.log_r);
  out.K = std::exp(par.log_K);
  out.q = std::exp(par.log_q);
  out.sigma_process = std::exp(par.log_sigma_process);
  out.sigma_index = std::exp(par.log_sigma_index);
  out.B0_frac = inv_logit(par.logit_B0_frac);
  out.B0 = out.B0_frac * out.K;
  return out;
}

double predicted_next_biomass(const double B, const double catch_value, const Transformed& tr) {
  const double production = tr.r * B * (1.0 - B / tr.K);
  return std::max(B + production - catch_value, 1e-9);
}

double joint_x(const ss::Data& data, const ss::Parameters& par, const Eigen::VectorXd& x) {
  const int n = static_cast<int>(data.catch_observed.size());
  const Transformed tr = transform(par);

  double nll = 0.0;
  nll += normal_nll(
      std::log(data.index_observed[0]) - (std::log(tr.q) + std::log(tr.B0)),
      tr.sigma_index);

  for (int t = 0; t < n - 1; ++t) {
    const double log_B_t = (t == 0) ? std::log(tr.B0) : x[t - 1];
    const double B_t = std::exp(log_B_t);
    const double pred_next =
        predicted_next_biomass(B_t, data.catch_observed[static_cast<std::size_t>(t)], tr);

    nll += normal_nll(x[t] - std::log(pred_next), tr.sigma_process);

    const double obs_residual =
        std::log(data.index_observed[static_cast<std::size_t>(t + 1)]) -
        (std::log(tr.q) + x[t]);

    nll += normal_nll(obs_residual, tr.sigma_index);
  }

  return nll;
}

Eigen::VectorXd deterministic_initial_x(const ss::Data& data, const ss::Parameters& par) {
  const int n = static_cast<int>(data.catch_observed.size());
  const Transformed tr = transform(par);
  Eigen::VectorXd x(n - 1);

  double B = tr.B0;
  for (int t = 0; t < n - 1; ++t) {
    B = predicted_next_biomass(B, data.catch_observed[static_cast<std::size_t>(t)], tr);
    x[t] = std::log(B);
  }
  return x;
}

Eigen::VectorXd fd_grad_x(const ss::Data& data, const ss::Parameters& par, const Eigen::VectorXd& x) {
  Eigen::VectorXd grad(x.size());
  for (int i = 0; i < x.size(); ++i) {
    const double h = 1e-5 * (1.0 + std::abs(x[i]));
    Eigen::VectorXd xp = x;
    Eigen::VectorXd xm = x;
    xp[i] += h;
    xm[i] -= h;
    grad[i] = (joint_x(data, par, xp) - joint_x(data, par, xm)) / (2.0 * h);
  }
  return grad;
}

class ObjX {
 public:
  ObjX(const ss::Data& data, const ss::Parameters& par) : data_(data), par_(par) {}
  double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad) {
    const double f = joint_x(data_, par_, x);
    grad = fd_grad_x(data_, par_, x);
    return f;
  }
 private:
  const ss::Data& data_;
  const ss::Parameters& par_;
};

Eigen::VectorXd optimize_x(const ss::Data& data, const ss::Parameters& par) {
  Eigen::VectorXd x = deterministic_initial_x(data, par);

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
  ObjX obj(data, par);
  double f = 0.0;

  try {
    solver.minimize(obj, x, f);
  } catch (...) {
    const double gnorm = fd_grad_x(data, par, x).norm();
    if (!(std::isfinite(joint_x(data, par, x)) && gnorm < 1e-3)) throw;
  }
  return x;
}

Eigen::SparseMatrix<double> fd_tridiagonal_hessian_xx(const ss::Data& data,
                                                      const ss::Parameters& par,
                                                      const Eigen::VectorXd& x) {
  const int n = static_cast<int>(x.size());
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(static_cast<std::size_t>(3 * n));

  for (int j = 0; j < n; ++j) {
    const double h = 1e-5 * (1.0 + std::abs(x[j]));
    Eigen::VectorXd xp = x;
    Eigen::VectorXd xm = x;
    xp[j] += h;
    xm[j] -= h;

    const Eigen::VectorXd gp = fd_grad_x(data, par, xp);
    const Eigen::VectorXd gm = fd_grad_x(data, par, xm);

    for (int i = std::max(0, j - 1); i <= std::min(n - 1, j + 1); ++i) {
      const double hij = (gp[i] - gm[i]) / (2.0 * h);
      if (std::abs(hij) > 1e-12) triplets.emplace_back(i, j, hij);
    }
  }

  Eigen::SparseMatrix<double> H(n, n);
  H.setFromTriplets(triplets.begin(), triplets.end());
  Eigen::SparseMatrix<double> Hsym = 0.5 * (H + Eigen::SparseMatrix<double>(H.transpose()));
  return Hsym;
}

double sparse_logdet_ldlt(const Eigen::SparseMatrix<double>& H) {
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
  ldlt.compute(H);
  if (ldlt.info() != Eigen::Success) throw std::runtime_error("Sparse LDLT failed");

  const auto& D = ldlt.vectorD();
  double logdet = 0.0;
  for (int i = 0; i < D.size(); ++i) {
    if (!(D[i] > 0.0)) throw std::runtime_error("Hxx not positive definite");
    logdet += std::log(D[i]);
  }
  return logdet;
}

struct EvalResult {
  double objective = 0.0;
  double joint = 0.0;
  double logdet = 0.0;
  double grad_norm = 0.0;
  int nnz = 0;
};

EvalResult eval(const ss::Data& data, const ss::Parameters& par) {
  EvalResult out;
  const Eigen::VectorXd xhat = optimize_x(data, par);
  out.joint = joint_x(data, par, xhat);
  out.grad_norm = fd_grad_x(data, par, xhat).norm();

  const Eigen::SparseMatrix<double> H = fd_tridiagonal_hessian_xx(data, par, xhat);
  out.nnz = static_cast<int>(H.nonZeros());
  out.logdet = sparse_logdet_ldlt(H);

  const double n_x = static_cast<double>(xhat.size());
  const double correction = 0.5 * out.logdet - 0.5 * n_x * std::log(2.0 * M_PI);
  out.objective = out.joint + correction;
  return out;
}

ss::Parameters make_par() {
  ss::Parameters par;
  par.log_r = std::log(0.5);
  par.log_K = std::log(700.0);
  par.log_q = std::log(0.0024);
  par.log_sigma_process = std::log(0.15);
  par.log_sigma_index = std::log(0.10);
  par.logit_B0_frac = std::log(0.90 / 0.10);
  return par;
}

}  // namespace

int main(int argc, char** argv) {
  int reps = 10;
  std::vector<int> lengths = {25, 50, 100, 250};

  if (argc > 1) reps = std::stoi(argv[1]);
  if (argc > 2) lengths = parse_lengths(argv[2]);

  const ss::Parameters par = make_par();

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "Quadra scaled latent-state tridiagonal Laplace benchmark\n";
  std::cout << "reps per n = " << reps << "\n\n";

  std::cout << std::setw(8) << "n"
            << std::setw(14) << "objective"
            << std::setw(14) << "joint"
            << std::setw(14) << "logdet"
            << std::setw(14) << "nnz"
            << std::setw(14) << "grad_norm"
            << std::setw(14) << "avg_ms"
            << "\n";

  for (const int n : lengths) {
    const ss::Data data = make_scaled_data(n);
    EvalResult last = eval(data, par);

    const auto t0 = Clock::now();
    for (int r = 0; r < reps; ++r) last = eval(data, par);
    const auto t1 = Clock::now();

    const double avg_ms = ms_between(t0, t1) / static_cast<double>(reps);

    std::cout << std::setw(8) << n
              << std::setw(14) << last.objective
              << std::setw(14) << last.joint
              << std::setw(14) << last.logdet
              << std::setw(14) << last.nnz
              << std::setw(14) << last.grad_norm
              << std::setw(14) << avg_ms
              << "\n";
  }

  return 0;
}
EOF

cat > run_quadra_scaled_latent_tridiagonal_benchmark.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-10}"
LENGTHS="${2:-25,50,100,250}"
CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/examples

set -x
"${CXX}" ${CXXFLAGS}   -Iexternal/Eigen   -Iexternal/LBFGSpp/include   -Iexamples/state_space_surplus_production   -Iexamples/surplus_production   examples/state_space_surplus_production/benchmark_latent_tridiagonal_scaled.cpp   -o build/examples/benchmark_latent_tridiagonal_scaled

./build/examples/benchmark_latent_tridiagonal_scaled "$REPS" "$LENGTHS"
EOF

chmod +x run_quadra_scaled_latent_tridiagonal_benchmark.sh

cat > examples/tmb_state_space_surplus/benchmark_scaled_fixed_theta_tmb.R <<'EOF'
#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[[1]]) else 10L
lengths <- if (length(args) >= 2) as.integer(strsplit(args[[2]], ",")[[1]]) else c(25L, 50L, 100L, 250L)

cat("TMB scaled fixed-theta state-space surplus benchmark\n")
cat("===================================================\n")
cat("reps per n =", reps, "\n\n")

if (!requireNamespace("TMB", quietly = TRUE)) {
  cat("TMB is not installed. Skipping.\n")
  quit(status = 0)
}

library(TMB)

template <- file.path("examples", "tmb_state_space_surplus", "state_space_surplus_tmb.cpp")
dynlib_name <- "state_space_surplus_tmb"

if (!file.exists(TMB::dynlib(file.path("examples", "tmb_state_space_surplus", dynlib_name)))) {
  TMB::compile(template, flags = "-O2")
}

dyn.load(TMB::dynlib(file.path("examples", "tmb_state_space_surplus", dynlib_name)))

make_scaled_data <- function(n) {
  r <- 0.5
  K <- 700
  q <- 0.0024
  B <- 0.90 * K
  catch_observed <- numeric(n)
  index_observed <- numeric(n)

  for (t0 in seq_len(n)) {
    t <- t0 - 1
    seasonal <- sin(2 * pi * t / 17)
    trend <- 1 + 0.10 * sin(2 * pi * t / 53)
    C <- 88 * trend + 18 * seasonal
    catch_observed[t0] <- max(40, C)

    obs_error <- 0.05 * sin(2 * pi * t / 11) +
      0.025 * cos(2 * pi * t / 7)

    index_observed[t0] <- q * B * exp(obs_error)

    if (t0 < n) {
      production <- r * B * (1 - B / K)
      B <- max(B + production - catch_observed[t0], 1e-9)
    }
  }

  list(catch_observed = catch_observed, index_observed = index_observed)
}

cat(sprintf("%8s%14s%14s\n", "n", "objective", "avg_ms"))

for (n in lengths) {
  data <- make_scaled_data(n)

  parameters <- list(
    log_r = log(0.5),
    log_K = log(700.0),
    log_q = log(0.0024),
    log_sigma_process = log(0.15),
    log_sigma_index = log(0.10),
    logit_B0_frac = log(0.90 / 0.10),
    u = rep(0, n - 1)
  )

  obj <- TMB::MakeADFun(
    data = data,
    parameters = parameters,
    random = "u",
    DLL = dynlib_name,
    silent = TRUE
  )

  last <- obj$fn()
  gc()

  t0 <- proc.time()
  for (i in seq_len(reps)) last <- obj$fn()
  t1 <- proc.time()

  avg_ms <- as.numeric((t1 - t0)[["elapsed"]]) * 1000 / reps

  cat(sprintf("%8d%14.6f%14.6f\n", n, last, avg_ms))
}
EOF

chmod +x examples/tmb_state_space_surplus/benchmark_scaled_fixed_theta_tmb.R

cat > run_tmb_scaled_fixed_theta_benchmark.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-10}"
LENGTHS="${2:-25,50,100,250}"
Rscript examples/tmb_state_space_surplus/benchmark_scaled_fixed_theta_tmb.R "$REPS" "$LENGTHS"
EOF

chmod +x run_tmb_scaled_fixed_theta_benchmark.sh

cat > run_quadra_vs_tmb_scaled_fixed_theta_benchmark.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-10}"
LENGTHS="${2:-25,50,100,250}"

echo "== Quadra scaled latent-state tridiagonal =="
./run_quadra_scaled_latent_tridiagonal_benchmark.sh "$REPS" "$LENGTHS"

echo
echo "== TMB scaled AD/Laplace =="
./run_tmb_scaled_fixed_theta_benchmark.sh "$REPS" "$LENGTHS"
EOF

chmod +x run_quadra_vs_tmb_scaled_fixed_theta_benchmark.sh

cat <<'EOF'

Installed scaled Quadra vs TMB benchmark.

Run:
  ./run_quadra_vs_tmb_scaled_fixed_theta_benchmark.sh 10 25,50,100,250

For a quicker smoke test:
  ./run_quadra_vs_tmb_scaled_fixed_theta_benchmark.sh 3 25,50

EOF
