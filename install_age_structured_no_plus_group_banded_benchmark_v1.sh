#!/usr/bin/env bash
set -euo pipefail

# install_age_structured_no_plus_group_banded_benchmark_v1.sh
#
# Adds a no-plus-group age-structured recruitment deviation benchmark.
#
# Why:
#   The previous banded age-structured benchmark did not match dense/TMB because
#   the plus group creates an effectively long-range dependency:
#
#     N_plus[t+1] += N_plus[t] * exp(-M)
#
#   Once a cohort enters the plus group, its effect persists indefinitely.
#   That makes H_xx dense-ish, not band-limited by n_ages.
#
# This benchmark removes the plus-group accumulator so each recruitment cohort
# exits after n_ages years. Then H_xx should be banded with bandwidth n_ages - 1.
#
# Adds:
#   - Quadra no-plus banded benchmark
#   - matching TMB no-plus benchmark
#   - combined runner
#
# Goal:
#   validate objective agreement, then inspect speed.

mkdir -p examples/age_structured_recruitment examples/tmb_age_structured_recruitment

cat > examples/age_structured_recruitment/benchmark_age_structured_no_plus_banded.cpp <<'EOF'
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <LBFGS.h>

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

double logistic(const double x) {
  return 1.0 / (1.0 + std::exp(-x));
}

double normal_nll(const double residual, const double sigma) {
  const double z = residual / sigma;
  return std::log(sigma) + 0.5 * std::log(2.0 * M_PI) + 0.5 * z * z;
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

struct Derived {
  double R0 = 0.0;
  double M = 0.0;
  double q = 0.0;
  double sigma_R = 0.0;
  double sigma_index = 0.0;
};

Derived transform(const Parameters& par) {
  Derived d;
  d.R0 = std::exp(par.log_R0);
  d.M = std::exp(par.log_M);
  d.q = std::exp(par.log_q);
  d.sigma_R = std::exp(par.log_sigma_R);
  d.sigma_index = std::exp(par.log_sigma_index);
  return d;
}

std::vector<double> selectivity(const int n_ages) {
  std::vector<double> sel(static_cast<std::size_t>(n_ages));
  for (int a = 0; a < n_ages; ++a) {
    sel[static_cast<std::size_t>(a)] =
        logistic((static_cast<double>(a + 1) - 4.0) / 0.8);
  }
  return sel;
}

std::vector<double> simulate_index(const int n_years,
                                   const int n_ages,
                                   const Parameters& par) {
  const Derived p = transform(par);
  const std::vector<double> sel = selectivity(n_ages);

  std::vector<double> N(static_cast<std::size_t>(n_ages), 0.0);
  for (int a = 0; a < n_ages; ++a) {
    N[static_cast<std::size_t>(a)] =
        p.R0 * std::exp(-p.M * static_cast<double>(a));
  }

  std::vector<double> index(static_cast<std::size_t>(n_years), 0.0);

  for (int y = 0; y < n_years; ++y) {
    double vulnerable = 0.0;
    for (int a = 0; a < n_ages; ++a) {
      vulnerable += sel[static_cast<std::size_t>(a)] *
                    N[static_cast<std::size_t>(a)];
    }

    const double obs_error =
        0.04 * std::sin(2.0 * M_PI * static_cast<double>(y) / 13.0) +
        0.02 * std::cos(2.0 * M_PI * static_cast<double>(y) / 7.0);

    index[static_cast<std::size_t>(y)] =
        p.q * vulnerable * std::exp(obs_error);

    std::vector<double> N_next(static_cast<std::size_t>(n_ages), 0.0);
    const double recruitment_dev =
        0.15 * std::sin(2.0 * M_PI * static_cast<double>(y) / 17.0);
    N_next[0] = p.R0 * std::exp(recruitment_dev);

    for (int a = 1; a < n_ages; ++a) {
      N_next[static_cast<std::size_t>(a)] =
          N[static_cast<std::size_t>(a - 1)] * std::exp(-p.M);
    }

    // No plus group: oldest cohort exits after terminal age.

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

double joint_x(const Data& data, const Parameters& par, const Eigen::VectorXd& x) {
  const Derived p = transform(par);
  const std::vector<double> sel = selectivity(data.n_ages);

  std::vector<double> N(static_cast<std::size_t>(data.n_ages), 0.0);
  for (int a = 0; a < data.n_ages; ++a) {
    N[static_cast<std::size_t>(a)] =
        p.R0 * std::exp(-p.M * static_cast<double>(a));
  }

  double nll = 0.0;

  for (int y = 0; y < data.n_years; ++y) {
    double vulnerable = 0.0;
    for (int a = 0; a < data.n_ages; ++a) {
      vulnerable += sel[static_cast<std::size_t>(a)] *
                    N[static_cast<std::size_t>(a)];
    }

    const double obs_residual =
        std::log(data.index_obs[static_cast<std::size_t>(y)]) -
        (std::log(p.q) + std::log(vulnerable));

    nll += normal_nll(obs_residual, p.sigma_index);
    nll += normal_nll(x[y], p.sigma_R);

    std::vector<double> N_next(static_cast<std::size_t>(data.n_ages), 0.0);
    N_next[0] = p.R0 * std::exp(x[y]);

    for (int a = 1; a < data.n_ages; ++a) {
      N_next[static_cast<std::size_t>(a)] =
          N[static_cast<std::size_t>(a - 1)] * std::exp(-p.M);
    }

    // No plus group.

    N.swap(N_next);
  }

  return nll;
}

Eigen::VectorXd fd_grad_x(const Data& data,
                          const Parameters& par,
                          const Eigen::VectorXd& x) {
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

class ObjX {
 public:
  ObjX(const Data& data, const Parameters& par)
      : data_(data), par_(par) {}

  double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad) {
    const double f = joint_x(data_, par_, x);
    grad = fd_grad_x(data_, par_, x);
    return f;
  }

 private:
  const Data& data_;
  const Parameters& par_;
};

Eigen::VectorXd optimize_x(const Data& data, const Parameters& par) {
  Eigen::VectorXd x = Eigen::VectorXd::Zero(data.n_years);

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
    if (!(std::isfinite(joint_x(data, par, x)) && gnorm < 1e-3)) {
      throw;
    }
  }

  return x;
}

Eigen::SparseMatrix<double> fd_banded_hessian_xx(const Data& data,
                                                 const Parameters& par,
                                                 const Eigen::VectorXd& x,
                                                 const int bandwidth) {
  const int n = static_cast<int>(x.size());
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(static_cast<std::size_t>(n * (2 * bandwidth + 1)));

  const double f0 = joint_x(data, par, x);

  for (int i = 0; i < n; ++i) {
    const double hi = 1e-4 * (1.0 + std::abs(x[i]));

    Eigen::VectorXd xp = x;
    Eigen::VectorXd xm = x;
    xp[i] += hi;
    xm[i] -= hi;

    const double hii = (joint_x(data, par, xp) - 2.0 * f0 +
                        joint_x(data, par, xm)) / (hi * hi);

    if (std::abs(hii) > 1e-12) {
      triplets.emplace_back(i, i, hii);
    }

    const int jmax = std::min(n - 1, i + bandwidth);
    for (int j = i + 1; j <= jmax; ++j) {
      const double hj = 1e-4 * (1.0 + std::abs(x[j]));

      Eigen::VectorXd xpp = x;
      Eigen::VectorXd xpm = x;
      Eigen::VectorXd xmp = x;
      Eigen::VectorXd xmm = x;

      xpp[i] += hi; xpp[j] += hj;
      xpm[i] += hi; xpm[j] -= hj;
      xmp[i] -= hi; xmp[j] += hj;
      xmm[i] -= hi; xmm[j] -= hj;

      const double hij =
          (joint_x(data, par, xpp) - joint_x(data, par, xpm) -
           joint_x(data, par, xmp) + joint_x(data, par, xmm)) /
          (4.0 * hi * hj);

      if (std::abs(hij) > 1e-12) {
        triplets.emplace_back(i, j, hij);
        triplets.emplace_back(j, i, hij);
      }
    }
  }

  Eigen::SparseMatrix<double> H(n, n);
  H.setFromTriplets(triplets.begin(), triplets.end());
  return H;
}

double sparse_logdet_ldlt(const Eigen::SparseMatrix<double>& H) {
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
  ldlt.compute(H);

  if (ldlt.info() != Eigen::Success) {
    throw std::runtime_error("sparse LDLT failed");
  }

  const auto& D = ldlt.vectorD();

  double logdet = 0.0;
  for (int i = 0; i < D.size(); ++i) {
    if (!(D[i] > 0.0)) {
      throw std::runtime_error("sparse Hxx is not positive definite");
    }
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

EvalResult eval(const Data& data,
                const Parameters& par,
                const int bandwidth) {
  EvalResult out;

  const Eigen::VectorXd xhat = optimize_x(data, par);
  out.joint = joint_x(data, par, xhat);
  out.grad_norm = fd_grad_x(data, par, xhat).norm();

  const Eigen::SparseMatrix<double> H = fd_banded_hessian_xx(data, par, xhat, bandwidth);
  out.nnz = static_cast<int>(H.nonZeros());
  out.logdet = sparse_logdet_ldlt(H);

  const double n_x = static_cast<double>(xhat.size());
  out.objective = out.joint + 0.5 * out.logdet -
                  0.5 * n_x * std::log(2.0 * M_PI);

  return out;
}

}  // namespace

int main(int argc, char** argv) {
  int reps = 10;
  std::vector<int> lengths = {25, 50, 100, 250, 500, 1000};
  int n_ages = 10;
  int bandwidth = 9;

  if (argc > 1) reps = std::stoi(argv[1]);
  if (argc > 2) lengths = parse_lengths(argv[2]);
  if (argc > 3) n_ages = std::stoi(argv[3]);
  if (argc > 4) bandwidth = std::stoi(argv[4]);

  const Parameters par;

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "Quadra no-plus age-structured banded Laplace benchmark\n";
  std::cout << "reps per n = " << reps
            << ", ages = " << n_ages
            << ", bandwidth = " << bandwidth
            << "\n\n";

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
    EvalResult last = eval(data, par, bandwidth);

    const auto t0 = Clock::now();
    for (int r = 0; r < reps; ++r) {
      last = eval(data, par, bandwidth);
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

cat > run_quadra_age_structured_no_plus_banded_benchmark.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-10}"
LENGTHS="${2:-25,50,100,250,500,1000}"
AGES="${3:-10}"
BANDWIDTH="${4:-9}"

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/examples

set -x
"${CXX}" ${CXXFLAGS} \
  -Iexternal/Eigen \
  -Iexternal/LBFGSpp/include \
  examples/age_structured_recruitment/benchmark_age_structured_no_plus_banded.cpp \
  -o build/examples/benchmark_age_structured_no_plus_banded

./build/examples/benchmark_age_structured_no_plus_banded "$REPS" "$LENGTHS" "$AGES" "$BANDWIDTH"
EOF

chmod +x run_quadra_age_structured_no_plus_banded_benchmark.sh

cat > examples/tmb_age_structured_recruitment/age_structured_recruitment_no_plus_tmb.cpp <<'EOF'
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
    N_next.setZero();

    N_next(0) = R0 * exp(x(y));

    for (int a = 1; a < n_ages; ++a) {
      N_next(a) = N(a - 1) * exp(-M);
    }

    // No plus group.

    N = N_next;
  }

  return nll;
}
EOF

cat > examples/tmb_age_structured_recruitment/benchmark_age_structured_no_plus_tmb.R <<'EOF'
#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
reps <- if (length(args) >= 1) as.integer(args[[1]]) else 10L
lengths <- if (length(args) >= 2) as.integer(strsplit(args[[2]], ",")[[1]]) else c(25L, 50L, 100L, 250L, 500L, 1000L)
n_ages <- if (length(args) >= 3) as.integer(args[[3]]) else 10L

cat("TMB no-plus age-structured recruitment benchmark\n")
cat("================================================\n")
cat("reps per n =", reps, ", ages =", n_ages, "\n\n")

if (!requireNamespace("TMB", quietly = TRUE)) {
  cat("TMB is not installed. Skipping.\n")
  quit(status = 0)
}

library(TMB)

template <- file.path("examples", "tmb_age_structured_recruitment", "age_structured_recruitment_no_plus_tmb.cpp")
dynlib_name <- "age_structured_recruitment_no_plus_tmb"

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

    # No plus group.

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

chmod +x examples/tmb_age_structured_recruitment/benchmark_age_structured_no_plus_tmb.R

cat > run_tmb_age_structured_no_plus_benchmark.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-10}"
LENGTHS="${2:-25,50,100,250,500,1000}"
AGES="${3:-10}"

Rscript examples/tmb_age_structured_recruitment/benchmark_age_structured_no_plus_tmb.R "$REPS" "$LENGTHS" "$AGES"
EOF

chmod +x run_tmb_age_structured_no_plus_benchmark.sh

cat > run_quadra_vs_tmb_age_structured_no_plus_benchmark.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-10}"
LENGTHS="${2:-25,50,100,250,500,1000}"
AGES="${3:-10}"
BANDWIDTH="${4:-9}"

echo "== Quadra no-plus age-structured banded Laplace =="
./run_quadra_age_structured_no_plus_banded_benchmark.sh "$REPS" "$LENGTHS" "$AGES" "$BANDWIDTH"

echo
echo "== TMB no-plus age-structured AD/Laplace =="
./run_tmb_age_structured_no_plus_benchmark.sh "$REPS" "$LENGTHS" "$AGES"
EOF

chmod +x run_quadra_vs_tmb_age_structured_no_plus_benchmark.sh

cat <<'EOF'

Installed no-plus age-structured banded benchmark.

Run:
  ./run_quadra_vs_tmb_age_structured_no_plus_benchmark.sh 10 25,50,100,250,500,1000 10 9

Quick:
  ./run_quadra_vs_tmb_age_structured_no_plus_benchmark.sh 3 25,50 10 9

Note:
  This removes the plus group so the recruitment-deviation Hessian is truly band-limited.
  The original plus-group model has long-range dependencies and needs a different sparse strategy.

EOF
