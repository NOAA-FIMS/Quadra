#!/usr/bin/env bash
set -euo pipefail

# install_age_structured_dense_laplace_validation_v1.sh
#
# Adds a correctness-first dense Laplace validation for the age-structured
# recruitment-deviation model.
#
# Why:
#   The previous Quadra benchmark used a diagonal Hessian approximation.
#   That is mathematically wrong for an age-structured model because a
#   recruitment deviation affects future numbers-at-age and future survey
#   observations. The H_xx structure is banded/dense-ish, not diagonal.
#
# This patch adds:
#   - full random-effects optimization with LBFGS++
#   - dense finite-difference H_xx
#   - dense Cholesky Laplace correction
#
# It is not meant to be fast. It is meant to match TMB first.
#
# Recommended validation:
#   ./run_quadra_age_structured_dense_laplace_validation.sh 3 25,50 10
#   ./run_tmb_age_structured_recruitment_benchmark.sh 3 25,50 10
#
# Once objectives match, the next step is a banded Hessian with bandwidth
# approximately n_ages.

mkdir -p examples/age_structured_recruitment

cat > examples/age_structured_recruitment/benchmark_age_structured_dense_laplace_validation.cpp <<'EOF'
#include <Eigen/Dense>
#include <LBFGS.h>

#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
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

    N_next[static_cast<std::size_t>(data.n_ages - 1)] +=
        N[static_cast<std::size_t>(data.n_ages - 1)] * std::exp(-p.M);

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

Eigen::MatrixXd fd_dense_hessian_xx(const Data& data,
                                    const Parameters& par,
                                    const Eigen::VectorXd& x) {
  const int n = static_cast<int>(x.size());
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(n, n);
  const double f0 = joint_x(data, par, x);

  for (int i = 0; i < n; ++i) {
    const double hi = 1e-4 * (1.0 + std::abs(x[i]));

    Eigen::VectorXd xp = x;
    Eigen::VectorXd xm = x;
    xp[i] += hi;
    xm[i] -= hi;

    H(i, i) = (joint_x(data, par, xp) - 2.0 * f0 +
               joint_x(data, par, xm)) / (hi * hi);

    for (int j = i + 1; j < n; ++j) {
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

      H(i, j) = hij;
      H(j, i) = hij;
    }
  }

  return H;
}

struct EvalResult {
  double objective = 0.0;
  double joint = 0.0;
  double logdet = 0.0;
  double grad_norm = 0.0;
  double min_eigenvalue = 0.0;
};

EvalResult eval(const Data& data, const Parameters& par) {
  EvalResult out;

  const Eigen::VectorXd xhat = optimize_x(data, par);
  out.joint = joint_x(data, par, xhat);
  out.grad_norm = fd_grad_x(data, par, xhat).norm();

  const Eigen::MatrixXd H = fd_dense_hessian_xx(data, par, xhat);

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H);
  out.min_eigenvalue = es.eigenvalues().minCoeff();

  Eigen::LLT<Eigen::MatrixXd> llt(H);
  if (llt.info() != Eigen::Success) {
    throw std::runtime_error("dense Hxx is not SPD");
  }

  const auto& L = llt.matrixL();
  out.logdet = 0.0;
  for (int i = 0; i < H.rows(); ++i) {
    out.logdet += 2.0 * std::log(L(i, i));
  }

  const double n_x = static_cast<double>(xhat.size());
  out.objective = out.joint + 0.5 * out.logdet -
                  0.5 * n_x * std::log(2.0 * M_PI);

  return out;
}

}  // namespace

int main(int argc, char** argv) {
  int reps = 3;
  std::vector<int> lengths = {25, 50};
  int n_ages = 10;

  if (argc > 1) reps = std::stoi(argv[1]);
  if (argc > 2) lengths = parse_lengths(argv[2]);
  if (argc > 3) n_ages = std::stoi(argv[3]);

  const Parameters par;

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "Quadra age-structured dense Laplace validation\n";
  std::cout << "reps per n = " << reps << ", ages = " << n_ages << "\n\n";

  std::cout << std::setw(8) << "n"
            << std::setw(14) << "objective"
            << std::setw(14) << "joint"
            << std::setw(14) << "logdet"
            << std::setw(14) << "grad_norm"
            << std::setw(14) << "min_eig"
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
              << std::setw(14) << last.grad_norm
              << std::setw(14) << last.min_eigenvalue
              << std::setw(14) << avg_ms
              << "\n";
  }

  return 0;
}
EOF

cat > run_quadra_age_structured_dense_laplace_validation.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-3}"
LENGTHS="${2:-25,50}"
AGES="${3:-10}"

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/examples

set -x
"${CXX}" ${CXXFLAGS} \
  -Iexternal/Eigen \
  -Iexternal/LBFGSpp/include \
  examples/age_structured_recruitment/benchmark_age_structured_dense_laplace_validation.cpp \
  -o build/examples/benchmark_age_structured_dense_laplace_validation

./build/examples/benchmark_age_structured_dense_laplace_validation "$REPS" "$LENGTHS" "$AGES"
EOF

chmod +x run_quadra_age_structured_dense_laplace_validation.sh

cat > run_validate_age_structured_against_tmb.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-3}"
LENGTHS="${2:-25,50}"
AGES="${3:-10}"

echo "== Quadra dense Laplace validation =="
./run_quadra_age_structured_dense_laplace_validation.sh "$REPS" "$LENGTHS" "$AGES"

echo
echo "== TMB AD/Laplace =="
./run_tmb_age_structured_recruitment_benchmark.sh "$REPS" "$LENGTHS" "$AGES"
EOF

chmod +x run_validate_age_structured_against_tmb.sh

cat <<'EOF'

Installed age-structured dense Laplace validation.

Run:
  ./run_validate_age_structured_against_tmb.sh 3 25,50 10

Do not use the old diagonal age benchmark for speed claims yet.
It is structurally approximate and should not match TMB.

EOF
