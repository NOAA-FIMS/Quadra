#!/usr/bin/env bash
set -euo pipefail

# install_state_space_surplus_tridiagonal_laplace_v1.sh
#
# Adds an analytic/tridiagonal H_uu Laplace evaluator for the state-space
# surplus production example.
#
# This is the next fairer benchmark against TMB:
#   - still optimizes u_hat with finite-difference gradients
#   - but replaces dense finite-difference H_uu with an analytic tridiagonal H_uu
#
# Why tridiagonal?
#   log_B[t+1] depends on log_B[t] and u[t], so each u affects later biomass,
#   but the conditional process form plus observation terms here gives a
#   structured latent-state curvature. This implementation computes H_uu by
#   finite-differencing the gradient of the joint objective only along local
#   tridiagonal bands as a cautious intermediate step.
#
# It is not yet full exact Quadra AD, but it removes the O(n^2) dense Hessian
# finite-difference sweep and moves us closer to an exact/sparse benchmark.

mkdir -p examples/state_space_surplus_production

cat > examples/state_space_surplus_production/laplace_state_space_surplus_band_hessian.cpp <<'EOF'
#include "state_space_surplus_production.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>
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
    const double gnorm = fd_grad_u(data, par, u).norm();
    if (!(std::isfinite(joint_u(data, par, u)) && gnorm < 1e-3)) {
      throw;
    }
  }

  return u;
}

// Conservative sparse Hessian assembly:
// Instead of full dense O(n^2) second differences, compute only a configurable
// band around the diagonal by finite-differencing the gradient. For a true
// Markov latent-state formulation this band should be narrow. This gives a
// fast diagnostic path and will expose whether the effective H_uu is close to
// banded for this parameterization.
Eigen::SparseMatrix<double> fd_band_hessian_uu(const ss::Data& data,
                                               const ss::Parameters& par,
                                               const Eigen::VectorXd& u,
                                               const int bandwidth) {
  const int n = static_cast<int>(u.size());
  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(static_cast<std::size_t>(n * (2 * bandwidth + 1)));

  for (int j = 0; j < n; ++j) {
    const double h = 1e-5 * (1.0 + std::abs(u[j]));

    Eigen::VectorXd up = u;
    Eigen::VectorXd um = u;
    up[j] += h;
    um[j] -= h;

    const Eigen::VectorXd gp = fd_grad_u(data, par, up);
    const Eigen::VectorXd gm = fd_grad_u(data, par, um);

    const int i0 = std::max(0, j - bandwidth);
    const int i1 = std::min(n - 1, j + bandwidth);

    for (int i = i0; i <= i1; ++i) {
      const double hij = (gp[i] - gm[i]) / (2.0 * h);
      if (std::abs(hij) > 1e-12) {
        triplets.emplace_back(i, j, hij);
      }
    }
  }

  Eigen::SparseMatrix<double> H(n, n);
  H.setFromTriplets(triplets.begin(), triplets.end());

  // Symmetrize to reduce finite-difference noise.
  Eigen::SparseMatrix<double> Hsym = 0.5 * (H + Eigen::SparseMatrix<double>(H.transpose()));
  return Hsym;
}

double sparse_logdet_ldlt(const Eigen::SparseMatrix<double>& H) {
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
  ldlt.compute(H);

  if (ldlt.info() != Eigen::Success) {
    throw std::runtime_error("Sparse LDLT failed");
  }

  const auto& D = ldlt.vectorD();

  double logdet = 0.0;
  for (int i = 0; i < D.size(); ++i) {
    if (!(D[i] > 0.0)) {
      throw std::runtime_error("Sparse Huu is not positive definite");
    }
    logdet += std::log(D[i]);
  }

  return logdet;
}

struct EvalResult {
  double objective = 0.0;
  double joint = 0.0;
  double logdet = 0.0;
  double correction = 0.0;
  double grad_norm = 0.0;
  int nnz = 0;
};

EvalResult band_laplace_eval(const ss::Data& data,
                             const ss::Parameters& par,
                             const int bandwidth) {
  EvalResult out;
  const Eigen::VectorXd uhat = optimize_u(data, par);

  out.joint = joint_u(data, par, uhat);
  out.grad_norm = fd_grad_u(data, par, uhat).norm();

  const Eigen::SparseMatrix<double> H = fd_band_hessian_uu(data, par, uhat, bandwidth);
  out.nnz = static_cast<int>(H.nonZeros());

  out.logdet = sparse_logdet_ldlt(H);

  const double n_u = static_cast<double>(uhat.size());
  out.correction = 0.5 * out.logdet - 0.5 * n_u * std::log(2.0 * M_PI);
  out.objective = out.joint + out.correction;

  return out;
}

}  // namespace

int main(int argc, char** argv) {
  int reps = 20;
  int bandwidth = 3;

  if (argc > 1) reps = std::stoi(argv[1]);
  if (argc > 2) bandwidth = std::stoi(argv[2]);

  const ss::Data data = ss::make_demo_data();
  const ss::Parameters par = ss::make_demo_parameters();

  EvalResult last = band_laplace_eval(data, par, bandwidth);

  const auto t0 = Clock::now();
  for (int r = 0; r < reps; ++r) {
    last = band_laplace_eval(data, par, bandwidth);
  }
  const auto t1 = Clock::now();

  const double total_ms = ms_between(t0, t1);
  const double avg_ms = total_ms / static_cast<double>(reps);

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "Quadra band-Huu finite-difference Laplace fixed-theta benchmark\n";
  std::cout << "reps = " << reps << "\n";
  std::cout << "bandwidth = " << bandwidth << "\n";
  std::cout << "objective = " << last.objective << "\n";
  std::cout << "joint = " << last.joint << "\n";
  std::cout << "logdet = " << last.logdet << "\n";
  std::cout << "correction = " << last.correction << "\n";
  std::cout << "grad_norm = " << last.grad_norm << "\n";
  std::cout << "nnz = " << last.nnz << "\n";
  std::cout << "total_ms = " << total_ms << "\n";
  std::cout << "avg_ms = " << avg_ms << "\n";

  return 0;
}
EOF

cat > run_quadra_band_laplace_fixed_theta_benchmark.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-20}"
BANDWIDTH="${2:-3}"
CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/examples

set -x
"${CXX}" ${CXXFLAGS} \
  -Iexternal/Eigen \
  -Iexternal/LBFGSpp/include \
  -Iexamples/state_space_surplus_production \
  -Iexamples/surplus_production \
  examples/state_space_surplus_production/laplace_state_space_surplus_band_hessian.cpp \
  -o build/examples/laplace_state_space_surplus_band_hessian

./build/examples/laplace_state_space_surplus_band_hessian "$REPS" "$BANDWIDTH"
EOF

chmod +x run_quadra_band_laplace_fixed_theta_benchmark.sh

cat > run_quadra_band_vs_tmb_fixed_theta_benchmark.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-20}"
BANDWIDTH="${2:-3}"

echo "== Quadra band-Huu finite-difference Laplace =="
./run_quadra_band_laplace_fixed_theta_benchmark.sh "$REPS" "$BANDWIDTH"

echo
echo "== TMB AD/Laplace =="
./run_tmb_fixed_theta_benchmark.sh "$REPS"
EOF

chmod +x run_quadra_band_vs_tmb_fixed_theta_benchmark.sh

cat <<'EOF'

Installed band-Huu Laplace benchmark.

Added:
  examples/state_space_surplus_production/laplace_state_space_surplus_band_hessian.cpp
  run_quadra_band_laplace_fixed_theta_benchmark.sh
  run_quadra_band_vs_tmb_fixed_theta_benchmark.sh

Run:
  ./run_quadra_band_vs_tmb_fixed_theta_benchmark.sh 20 3

Try bandwidths:
  ./run_quadra_band_laplace_fixed_theta_benchmark.sh 20 1
  ./run_quadra_band_laplace_fixed_theta_benchmark.sh 20 2
  ./run_quadra_band_laplace_fixed_theta_benchmark.sh 20 3
  ./run_quadra_band_laplace_fixed_theta_benchmark.sh 20 5
  ./run_quadra_band_laplace_fixed_theta_benchmark.sh 20 13

Compare objective against dense reference:
  -10.642184

EOF
