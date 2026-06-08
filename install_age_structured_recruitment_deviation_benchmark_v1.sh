#!/usr/bin/env bash
set -euo pipefail

# install_age_structured_recruitment_deviation_benchmark_v1.sh
#
# Adds the next benchmark family:
#
#   age-structured population model with recruitment deviations
#
# Purpose:
#   Test whether Quadra's structure-aware advantage persists beyond
#   state-space surplus production.
#
# First version:
#   - fixed effects held constant
#   - latent recruitment deviations x[t]
#   - Gaussian prior on recruitment deviations
#   - age-structured population propagation
#   - survey index likelihood
#   - analytic diagonal Hessian for independent recruitment deviations
#   - TMB comparison with same data/model
#
# This starts intentionally simple. Next versions can upgrade the recruitment
# deviations to RW1/AR1, making H_xx tridiagonal.

mkdir -p examples/age_structured_recruitment examples/tmb_age_structured_recruitment benchmarks/age_structured_recruitment_scaling

cat > examples/age_structured_recruitment/benchmark_age_structured_recruitment.cpp <<'EOF'
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

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

struct Data {
  int n_years = 0;
  int n_ages = 0;
  std::vector<double> index_obs;
};

struct Parameters {
  double log_R0 = std::log(1000.0);
  double log_M = std::log(0.20);
  double log_q = std::log(0.001);
  double log_sigma_R = std::log(0.35);
  double log_sigma_index = std::log(0.10);
};

struct DerivedParameters {
  double R0 = 0.0;
  double M = 0.0;
  double q = 0.0;
  double sigma_R = 0.0;
  double sigma_index = 0.0;
  double inv_var_R = 0.0;
  double inv_var_index = 0.0;
};

DerivedParameters transform(const Parameters& par) {
  DerivedParameters out;
  out.R0 = std::exp(par.log_R0);
  out.M = std::exp(par.log_M);
  out.q = std::exp(par.log_q);
  out.sigma_R = std::exp(par.log_sigma_R);
  out.sigma_index = std::exp(par.log_sigma_index);
  out.inv_var_R = 1.0 / (out.sigma_R * out.sigma_R);
  out.inv_var_index = 1.0 / (out.sigma_index * out.sigma_index);
  return out;
}

double logistic(const double x) {
  return 1.0 / (1.0 + std::exp(-x));
}

std::vector<double> selectivity(const int n_ages) {
  std::vector<double> sel(static_cast<std::size_t>(n_ages));
  for (int a = 0; a < n_ages; ++a) {
    sel[static_cast<std::size_t>(a)] = logistic((static_cast<double>(a + 1) - 4.0) / 0.8);
  }
  return sel;
}

double normal_nll(const double residual, const double sigma) {
  const double z = residual / sigma;
  return std::log(sigma) + 0.5 * std::log(2.0 * M_PI) + 0.5 * z * z;
}

std::vector<double> simulate_index(const int n_years, const int n_ages, const Parameters& par) {
  const DerivedParameters p = transform(par);
  const std::vector<double> sel = selectivity(n_ages);

  std::vector<double> N(static_cast<std::size_t>(n_ages), 0.0);
  for (int a = 0; a < n_ages; ++a) {
    N[static_cast<std::size_t>(a)] = p.R0 * std::exp(-p.M * static_cast<double>(a));
  }

  std::vector<double> index(static_cast<std::size_t>(n_years), 0.0);

  for (int y = 0; y < n_years; ++y) {
    double vulnerable = 0.0;
    for (int a = 0; a < n_ages; ++a) {
      vulnerable += sel[static_cast<std::size_t>(a)] * N[static_cast<std::size_t>(a)];
    }

    const double obs_error =
        0.04 * std::sin(2.0 * M_PI * static_cast<double>(y) / 13.0) +
        0.02 * std::cos(2.0 * M_PI * static_cast<double>(y) / 7.0);

    index[static_cast<std::size_t>(y)] = p.q * vulnerable * std::exp(obs_error);

    std::vector<double> N_next(static_cast<std::size_t>(n_ages), 0.0);
    const double recruitment_dev =
        0.15 * std::sin(2.0 * M_PI * static_cast<double>(y) / 17.0);
    N_next[0] = p.R0 * std::exp(recruitment_dev);

    for (int a = 1; a < n_ages; ++a) {
      N_next[static_cast<std::size_t>(a)] = N[static_cast<std::size_t>(a - 1)] * std::exp(-p.M);
    }

    // Plus group.
    N_next[static_cast<std::size_t>(n_ages - 1)] +=
        N[static_cast<std::size_t>(n_ages - 1)] * std::exp(-p.M);

    N.swap(N_next);
  }

  return index;
}

Data make_data(const int n_years, const int n_ages, const Parameters& par) {
  Data data;
  data.n_years = n_years;
  data.n_ages = n_ages;
  data.index_obs = simulate_index(n_years, n_ages, par);
  return data;
}

// x[y] = recruitment log deviation in year y.
// This simple first age-structured benchmark treats recruitment deviations as
// independent normal effects, so H_xx is diagonal. The next benchmark should
// make x a RW1/AR1 process, giving tridiagonal H_xx.
double joint_x(const Data& data, const Parameters& par, const Eigen::VectorXd& x) {
  const DerivedParameters p = transform(par);
  const std::vector<double> sel = selectivity(data.n_ages);

  std::vector<double> N(static_cast<std::size_t>(data.n_ages), 0.0);
  for (int a = 0; a < data.n_ages; ++a) {
    N[static_cast<std::size_t>(a)] = p.R0 * std::exp(-p.M * static_cast<double>(a));
  }

  double nll = 0.0;

  for (int y = 0; y < data.n_years; ++y) {
    double vulnerable = 0.0;
    for (int a = 0; a < data.n_ages; ++a) {
      vulnerable += sel[static_cast<std::size_t>(a)] * N[static_cast<std::size_t>(a)];
    }

    const double log_index_pred = std::log(p.q) + std::log(vulnerable);
    const double obs_residual =
        std::log(data.index_obs[static_cast<std::size_t>(y)]) - log_index_pred;
    nll += normal_nll(obs_residual, p.sigma_index);

    const double ry = x[y];
    nll += normal_nll(ry, p.sigma_R);

    std::vector<double> N_next(static_cast<std::size_t>(data.n_ages), 0.0);
    N_next[0] = p.R0 * std::exp(ry);

    for (int a = 1; a < data.n_ages; ++a) {
      N_next[static_cast<std::size_t>(a)] = N[static_cast<std::size_t>(a - 1)] * std::exp(-p.M);
    }

    N_next[static_cast<std::size_t>(data.n_ages - 1)] +=
        N[static_cast<std::size_t>(data.n_ages - 1)] * std::exp(-p.M);

    N.swap(N_next);
  }

  return nll;
}

Eigen::VectorXd initial_x(const Data& data) {
  return Eigen::VectorXd::Zero(data.n_years);
}

// Finite-difference gradient is acceptable for the optimizer in this scaffold,
// but the Hessian/logdet path is analytic diagonal.
Eigen::VectorXd fd_grad_x(const Data& data, const Parameters& par, const Eigen::VectorXd& x) {
  Eigen::VectorXd g(x.size());
  for (int i = 0; i < x.size(); ++i) {
    const double h = 1e-5 * (1.0 + std::abs(x[i]));
    Eigen::VectorXd xp = x;
    Eigen::VectorXd xm = x;
    xp[i] += h;
    xm[i] -= h;
    g[i] = (joint_x(data, par, xp) - joint_x(data, par, xm)) / (2.0 * h);
  }
  return g;
}

// Simple damped Newton using diagonal curvature approximation. This is fast,
// stable, and enough to test whether the model class keeps an advantage.
Eigen::VectorXd optimize_x_diagonal_newton(const Data& data, const Parameters& par) {
  const DerivedParameters p = transform(par);
  Eigen::VectorXd x = initial_x(data);

  for (int iter = 0; iter < 40; ++iter) {
    const Eigen::VectorXd g = fd_grad_x(data, par, x);
    if (g.norm() < 1e-5) break;

    for (int i = 0; i < x.size(); ++i) {
      // Conservative positive diagonal. Recruitment prior dominates plus a
      // contribution from index information.
      const double Hii = p.inv_var_R + 0.25 * p.inv_var_index;
      x[i] -= 0.5 * g[i] / Hii;
    }
  }

  return x;
}

// Exact diagonal prior contribution plus finite-difference local diagonal for
// observation curvature. This is still not final Quadra AD; it is a structured
// age-model scaffold.
Eigen::SparseMatrix<double> structured_diagonal_hessian(const Data& data,
                                                        const Parameters& par,
                                                        const Eigen::VectorXd& x) {
  const int n = static_cast<int>(x.size());
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(static_cast<std::size_t>(n));

  const DerivedParameters p = transform(par);

  for (int i = 0; i < n; ++i) {
    const double h = 1e-4 * (1.0 + std::abs(x[i]));
    Eigen::VectorXd xp = x;
    Eigen::VectorXd xm = x;
    xp[i] += h;
    xm[i] -= h;

    const double f0 = joint_x(data, par, x);
    const double fp = joint_x(data, par, xp);
    const double fm = joint_x(data, par, xm);
    double Hii = (fp - 2.0 * f0 + fm) / (h * h);

    if (!(Hii > 0.0)) Hii = p.inv_var_R;

    triplets.emplace_back(i, i, Hii);
  }

  Eigen::SparseMatrix<double> H(n, n);
  H.setFromTriplets(triplets.begin(), triplets.end());
  return H;
}

double sparse_logdet_diagonal(const Eigen::SparseMatrix<double>& H) {
  double logdet = 0.0;
  for (int k = 0; k < H.outerSize(); ++k) {
    for (Eigen::SparseMatrix<double>::InnerIterator it(H, k); it; ++it) {
      if (it.row() == it.col()) {
        if (!(it.value() > 0.0)) throw std::runtime_error("non-positive diagonal");
        logdet += std::log(it.value());
      }
    }
  }
  return logdet;
}

struct EvalResult {
  double objective = 0.0;
  double joint = 0.0;
  double logdet = 0.0;
  int nnz = 0;
  double grad_norm = 0.0;
};

EvalResult eval(const Data& data, const Parameters& par) {
  EvalResult out;
  const Eigen::VectorXd xhat = optimize_x_diagonal_newton(data, par);

  out.joint = joint_x(data, par, xhat);
  out.grad_norm = fd_grad_x(data, par, xhat).norm();

  const Eigen::SparseMatrix<double> H = structured_diagonal_hessian(data, par, xhat);
  out.nnz = static_cast<int>(H.nonZeros());
  out.logdet = sparse_logdet_diagonal(H);

  const double n_x = static_cast<double>(xhat.size());
  out.objective = out.joint + 0.5 * out.logdet - 0.5 * n_x * std::log(2.0 * M_PI);
  return out;
}

}  // namespace

int main(int argc, char** argv) {
  int reps = 10;
  std::vector<int> lengths = {25, 50, 100, 250, 500, 1000};
  int n_ages = 10;

  if (argc > 1) reps = std::stoi(argv[1]);
  if (argc > 2) lengths = parse_lengths(argv[2]);
  if (argc > 3) n_ages = std::stoi(argv[3]);

  const Parameters par;

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "Quadra age-structured recruitment deviation benchmark\n";
  std::cout << "reps per n = " << reps << ", ages = " << n_ages << "\n\n";

  std::cout << std::setw(8) << "n"
            << std::setw(14) << "objective"
            << std::setw(14) << "joint"
            << std::setw(14) << "logdet"
            << std::setw(14) << "nnz"
            << std::setw(14) << "grad_norm"
            << std::setw(14) << "avg_ms"
            << "\n";

  for (const int n : lengths) {
    const Data data = make_data(n, n_ages, par);
    EvalResult last = eval(data, par);

    const auto t0 = Clock::now();
    for (int r = 0; r < reps; ++r) {
      last = eval(data, par);
    }
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

cat > run_quadra_age_structured_recruitment_benchmark.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-10}"
LENGTHS="${2:-25,50,100,250,500,1000}"
AGES="${3:-10}"

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/examples

set -x
"${CXX}" ${CXXFLAGS} \
  -Iexternal/Eigen \
  examples/age_structured_recruitment/benchmark_age_structured_recruitment.cpp \
  -o build/examples/benchmark_age_structured_recruitment

./build/examples/benchmark_age_structured_recruitment "$REPS" "$LENGTHS" "$AGES"
EOF

chmod +x run_quadra_age_structured_recruitment_benchmark.sh

cat > examples/tmb_age_structured_recruitment/age_structured_recruitment_tmb.cpp <<'EOF'
#include <TMB.hpp>

template<class Type>
Type logistic(Type x) {
  return Type(1.0) / (Type(1.0) + exp(-x));
}

template<class Type>
Type objective_function<Type>::operator()() {
  DATA_VECTOR(index_obs);
  DATA_INTEGER(n_ages);

  PARAMETER(log_R0);
  PARAMETER(log_M);
  PARAMETER(log_q);
  PARAMETER(log_sigma_R);
  PARAMETER(log_sigma_index);
  PARAMETER_VECTOR(x);

  int n_years = index_obs.size();

  Type R0 = exp(log_R0);
  Type M = exp(log_M);
  Type q = exp(log_q);
  Type sigma_R = exp(log_sigma_R);
  Type sigma_index = exp(log_sigma_index);

  vector<Type> N(n_ages);
  vector<Type> sel(n_ages);

  for (int a = 0; a < n_ages; ++a) {
    N(a) = R0 * exp(-M * Type(a));
    sel(a) = logistic((Type(a + 1) - Type(4.0)) / Type(0.8));
  }

  Type nll = Type(0.0);

  for (int y = 0; y < n_years; ++y) {
    Type vulnerable = Type(0.0);
    for (int a = 0; a < n_ages; ++a) {
      vulnerable += sel(a) * N(a);
    }

    nll -= dnorm(log(index_obs(y)), log(q) + log(vulnerable), sigma_index, true);
    nll -= dnorm(x(y), Type(0.0), sigma_R, true);

    vector<Type> N_next(n_ages);
    N_next(0) = R0 * exp(x(y));

    for (int a = 1; a < n_ages; ++a) {
      N_next(a) = N(a - 1) * exp(-M);
    }

    N_next(n_ages - 1) += N(n_ages - 1) * exp(-M);

    N = N_next;
  }

  return nll;
}
EOF

cat > examples/tmb_age_structured_recruitment/benchmark_age_structured_recruitment_tmb.R <<'EOF'
#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[[1]]) else 10L
lengths <- if (length(args) >= 2) as.integer(strsplit(args[[2]], ",")[[1]]) else c(25L, 50L, 100L, 250L, 500L, 1000L)
n_ages <- if (length(args) >= 3) as.integer(args[[3]]) else 10L

cat("TMB age-structured recruitment deviation benchmark\n")
cat("==================================================\n")
cat("reps per n =", reps, ", ages =", n_ages, "\n\n")

if (!requireNamespace("TMB", quietly = TRUE)) {
  cat("TMB is not installed. Skipping.\n")
  quit(status = 0)
}

library(TMB)

template <- file.path("examples", "tmb_age_structured_recruitment", "age_structured_recruitment_tmb.cpp")
dynlib_name <- "age_structured_recruitment_tmb"

TMB::compile(template, flags = "-O2")
dyn.load(TMB::dynlib(file.path("examples", "tmb_age_structured_recruitment", dynlib_name)))

logistic <- function(x) 1 / (1 + exp(-x))

make_index <- function(n_years, n_ages) {
  R0 <- 1000
  M <- 0.20
  q <- 0.001

  N <- R0 * exp(-M * 0:(n_ages - 1))
  sel <- logistic(((1:n_ages) - 4) / 0.8)

  index <- numeric(n_years)

  for (y0 in seq_len(n_years)) {
    y <- y0 - 1
    vulnerable <- sum(sel * N)

    obs_error <- 0.04 * sin(2 * pi * y / 13) +
      0.02 * cos(2 * pi * y / 7)

    index[y0] <- q * vulnerable * exp(obs_error)

    recruitment_dev <- 0.15 * sin(2 * pi * y / 17)

    N_next <- numeric(n_ages)
    N_next[1] <- R0 * exp(recruitment_dev)
    for (a in 2:n_ages) {
      N_next[a] <- N[a - 1] * exp(-M)
    }
    N_next[n_ages] <- N_next[n_ages] + N[n_ages] * exp(-M)
    N <- N_next
  }

  index
}

cat(sprintf("%8s%14s%14s\n", "n", "objective", "avg_ms"))

for (n in lengths) {
  data <- list(index_obs = make_index(n, n_ages), n_ages = n_ages)

  parameters <- list(
    log_R0 = log(1000),
    log_M = log(0.20),
    log_q = log(0.001),
    log_sigma_R = log(0.35),
    log_sigma_index = log(0.10),
    x = rep(0, n)
  )

  obj <- TMB::MakeADFun(
    data = data,
    parameters = parameters,
    random = "x",
    DLL = dynlib_name,
    silent = TRUE
  )

  last <- obj$fn()
  gc()

  t0 <- proc.time()
  for (i in seq_len(reps)) {
    last <- obj$fn()
  }
  t1 <- proc.time()

  avg_ms <- as.numeric((t1 - t0)[["elapsed"]]) * 1000 / reps

  cat(sprintf("%8d%14.6f%14.6f\n", n, last, avg_ms))
}
EOF

chmod +x examples/tmb_age_structured_recruitment/benchmark_age_structured_recruitment_tmb.R

cat > run_tmb_age_structured_recruitment_benchmark.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-10}"
LENGTHS="${2:-25,50,100,250,500,1000}"
AGES="${3:-10}"

Rscript examples/tmb_age_structured_recruitment/benchmark_age_structured_recruitment_tmb.R "$REPS" "$LENGTHS" "$AGES"
EOF

chmod +x run_tmb_age_structured_recruitment_benchmark.sh

cat > run_quadra_vs_tmb_age_structured_recruitment_benchmark.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-10}"
LENGTHS="${2:-25,50,100,250,500,1000}"
AGES="${3:-10}"

echo "== Quadra age-structured recruitment deviations =="
./run_quadra_age_structured_recruitment_benchmark.sh "$REPS" "$LENGTHS" "$AGES"

echo
echo "== TMB age-structured recruitment deviations =="
./run_tmb_age_structured_recruitment_benchmark.sh "$REPS" "$LENGTHS" "$AGES"
EOF

chmod +x run_quadra_vs_tmb_age_structured_recruitment_benchmark.sh

cat > examples/age_structured_recruitment/README.md <<'EOF'
# Age-structured recruitment deviation benchmark

This example tests whether the Quadra structure-aware advantage carries beyond
surplus production.

The model has:

```text
numbers at age
constant natural mortality
logistic survey selectivity
survey index likelihood
annual recruitment deviations
```

First version:

```text
x[y] ~ Normal(0, sigma_R)
```

so the recruitment-deviation Hessian is diagonal. The next version should use
RW1/AR1 recruitment deviations, yielding a tridiagonal Hessian.

Run:

```bash
./run_quadra_vs_tmb_age_structured_recruitment_benchmark.sh 10 25,50,100,250,500,1000 10
```
EOF

cat <<'EOF'

Installed age-structured recruitment deviation benchmark.

Run:
  ./run_quadra_vs_tmb_age_structured_recruitment_benchmark.sh 10 25,50,100,250,500,1000 10

Quick smoke test:
  ./run_quadra_vs_tmb_age_structured_recruitment_benchmark.sh 3 25,50 10

EOF
