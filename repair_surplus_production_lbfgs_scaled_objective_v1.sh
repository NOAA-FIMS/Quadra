#!/usr/bin/env bash
set -euo pipefail

# repair_surplus_production_lbfgs_scaled_objective_v1.sh
#
# Fixes LBFGS++ line-search failure in the surplus production example by:
#   - optimizing centered/scaled variables x instead of raw log parameters
#   - using a smooth objective with soft priors instead of hard-ish guard penalties
#   - using a slightly larger finite-difference step for stable numerical gradients
#   - starting near the biologically reasonable region
#
# This still uses finite-difference gradients. The next step is Quadra AD.

mkdir -p .quadra_patch_backups

target="examples/surplus_production/fit_surplus_production_lbfgs.cpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  echo "Run add_missing_surplus_production_lbfgs_example_v1.sh first."
  exit 1
fi

cp "$target" ".quadra_patch_backups/fit_surplus_production_lbfgs.cpp.scaled_objective.$(date +%Y%m%d_%H%M%S).bak"

cat > "$target" <<'EOF'
#include "surplus_production.hpp"

#include <Eigen/Dense>
#include <LBFGS.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace sp = quadra_examples::surplus_production;

namespace {

struct TransformCenter {
  double log_r = std::log(0.35);
  double log_K = std::log(1250.0);
  double log_q = std::log(0.00145);
  double log_sigma_index = std::log(0.15);
  double logit_B0_frac = std::log(0.85 / 0.15);
};

sp::Parameters unpack_scaled(const Eigen::VectorXd& x,
                             const TransformCenter& center) {
  sp::Parameters par;
  par.log_r = center.log_r + x[0];
  par.log_K = center.log_K + x[1];
  par.log_q = center.log_q + x[2];
  par.log_sigma_index = center.log_sigma_index + x[3];
  par.logit_B0_frac = center.logit_B0_frac + x[4];
  return par;
}

Eigen::VectorXd pack_scaled(const sp::Parameters& par,
                            const TransformCenter& center) {
  Eigen::VectorXd x(5);
  x[0] = par.log_r - center.log_r;
  x[1] = par.log_K - center.log_K;
  x[2] = par.log_q - center.log_q;
  x[3] = par.log_sigma_index - center.log_sigma_index;
  x[4] = par.logit_B0_frac - center.logit_B0_frac;
  return x;
}

double square(const double x) { return x * x; }

double soft_prior_penalty(const sp::Derived& d,
                          const sp::Parameters& par,
                          const TransformCenter& center) {
  // These are weak regularization terms for a tiny toy dataset.
  // They stabilize K/q/r tradeoffs without forcing a single answer.
  double p = 0.0;

  p += 0.02 * square(par.log_r - center.log_r);
  p += 0.02 * square(par.log_K - center.log_K);
  p += 0.02 * square(par.log_q - center.log_q);
  p += 0.02 * square(par.logit_B0_frac - center.logit_B0_frac);

  // Prevent sigma from collapsing to zero in the toy example.
  const double log_sigma_floor = std::log(0.06);
  if (par.log_sigma_index < log_sigma_floor) {
    p += 20.0 * square(par.log_sigma_index - log_sigma_floor);
  }

  // Keep the deterministic trajectory away from numerical collapse.
  if (d.depletion_terminal < 0.05) {
    p += 100.0 * square(0.05 - d.depletion_terminal);
  }

  return p;
}

double objective_value(const sp::Data& data,
                       const TransformCenter& center,
                       const Eigen::VectorXd& x) {
  try {
    const sp::Parameters par = unpack_scaled(x, center);
    const sp::Derived d = sp::evaluate_derived(data, par);

    if (!std::isfinite(d.r) || !std::isfinite(d.K) ||
        !std::isfinite(d.q) || !std::isfinite(d.sigma_index)) {
      return std::numeric_limits<double>::infinity();
    }

    return sp::negative_log_likelihood(data, par) +
           soft_prior_penalty(d, par, center);
  } catch (...) {
    return std::numeric_limits<double>::infinity();
  }
}

Eigen::VectorXd finite_difference_gradient(const sp::Data& data,
                                           const TransformCenter& center,
                                           const Eigen::VectorXd& x) {
  Eigen::VectorXd grad(x.size());

  for (int i = 0; i < x.size(); ++i) {
    const double step = 1e-4 * (1.0 + std::abs(x[i]));

    Eigen::VectorXd plus = x;
    Eigen::VectorXd minus = x;
    plus[i] += step;
    minus[i] -= step;

    const double f_plus = objective_value(data, center, plus);
    const double f_minus = objective_value(data, center, minus);

    if (!std::isfinite(f_plus) || !std::isfinite(f_minus)) {
      grad[i] = 0.0;
    } else {
      grad[i] = (f_plus - f_minus) / (2.0 * step);
    }
  }

  return grad;
}

class Objective {
 public:
  Objective(const sp::Data& data, const TransformCenter& center)
      : data_(data), center_(center) {}

  double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad) {
    const double f = objective_value(data_, center_, x);
    grad = finite_difference_gradient(data_, center_, x);
    return f;
  }

 private:
  const sp::Data& data_;
  const TransformCenter& center_;
};

void print_parameter_comparison(const sp::Data& data,
                                const sp::Parameters& initial,
                                const sp::Parameters& estimated) {
  const sp::Derived d0 = sp::evaluate_derived(data, initial);
  const sp::Derived d1 = sp::evaluate_derived(data, estimated);

  std::cout << "\nParameter comparison\n";
  std::cout << std::setw(18) << "quantity"
            << std::setw(16) << "initial"
            << std::setw(16) << "estimated"
            << "\n";

  std::cout << std::setw(18) << "r"
            << std::setw(16) << d0.r
            << std::setw(16) << d1.r
            << "\n";

  std::cout << std::setw(18) << "K"
            << std::setw(16) << d0.K
            << std::setw(16) << d1.K
            << "\n";

  std::cout << std::setw(18) << "q"
            << std::setw(16) << d0.q
            << std::setw(16) << d1.q
            << "\n";

  std::cout << std::setw(18) << "sigma_index"
            << std::setw(16) << d0.sigma_index
            << std::setw(16) << d1.sigma_index
            << "\n";

  std::cout << std::setw(18) << "B0/K"
            << std::setw(16) << d0.B0_frac
            << std::setw(16) << d1.B0_frac
            << "\n";

  std::cout << std::setw(18) << "MSY"
            << std::setw(16) << d0.MSY
            << std::setw(16) << d1.MSY
            << "\n";

  std::cout << std::setw(18) << "B_terminal/K"
            << std::setw(16) << d0.depletion_terminal
            << std::setw(16) << d1.depletion_terminal
            << "\n";
}

}  // namespace

int main() {
  const sp::Data data = sp::make_demo_data();
  const TransformCenter center;

  sp::Parameters initial;
  initial.log_r = std::log(0.28);
  initial.log_K = std::log(1800.0);
  initial.log_q = std::log(0.0010);
  initial.log_sigma_index = std::log(0.25);
  initial.logit_B0_frac = std::log(0.80 / 0.20);

  Eigen::VectorXd x = pack_scaled(initial, center);

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "Fit surplus production example with LBFGS++\n";
  std::cout << "==========================================\n\n";
  std::cout << "initial objective = "
            << objective_value(data, center, x) << "\n\n";

  LBFGSpp::LBFGSParam<double> param;
  param.epsilon = 1e-5;
  param.max_iterations = 500;
  param.max_linesearch = 100;
  param.m = 8;
  param.ftol = 1e-4;
  param.wolfe = 0.9;
  param.min_step = 1e-20;
  param.max_step = 1.0;

  LBFGSpp::LBFGSSolver<double> solver(param);
  Objective objective(data, center);

  double final_objective = 0.0;
  int iterations = 0;

  try {
    iterations = solver.minimize(objective, x, final_objective);
  } catch (const std::exception& e) {
    std::cerr << "LBFGS++ failed: " << e.what() << "\n";
    std::cerr << "Current objective = "
              << objective_value(data, center, x) << "\n";
    std::cerr << "Current gradient norm = "
              << finite_difference_gradient(data, center, x).norm() << "\n";
    return 1;
  }

  const sp::Parameters estimated = unpack_scaled(x, center);
  const Eigen::VectorXd final_grad = finite_difference_gradient(data, center, x);

  std::cout << "Fit summary\n";
  std::cout << "  iterations = " << iterations << "\n";
  std::cout << "  objective  = " << final_objective << "\n";
  std::cout << "  grad_norm  = " << final_grad.norm() << "\n";

  print_parameter_comparison(data, initial, estimated);

  std::cout << "\nEstimated model report\n";
  std::cout << "======================\n";
  sp::print_report(data, estimated);

  return 0;
}
EOF

cat <<'EOF'

Patched LBFGS++ surplus production example with scaled variables.

Run:
  ./run_fit_surplus_production_lbfgs_example.sh

EOF
