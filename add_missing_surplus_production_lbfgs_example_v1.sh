#!/usr/bin/env bash
set -euo pipefail

# add_missing_surplus_production_lbfgs_example_v1.sh
#
# Adds the missing LBFGS++ surplus production fitting example.
#
# Existing files expected:
#   examples/surplus_production/surplus_production.hpp
#
# Adds:
#   examples/surplus_production/fit_surplus_production_lbfgs.cpp
#   run_fit_surplus_production_lbfgs_example.sh

mkdir -p examples/surplus_production

if [[ ! -f examples/surplus_production/surplus_production.hpp ]]; then
  echo "ERROR: missing examples/surplus_production/surplus_production.hpp"
  exit 1
fi

if [[ ! -d external/Eigen ]]; then
  echo "ERROR: missing external/Eigen"
  exit 1
fi

if [[ ! -d external/LBFGSpp/include ]]; then
  echo "ERROR: missing external/LBFGSpp/include"
  exit 1
fi

cat > examples/surplus_production/fit_surplus_production_lbfgs.cpp <<'EOF'
#include "surplus_production.hpp"

#include <Eigen/Dense>
#include <LBFGS.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

namespace sp = quadra_examples::surplus_production;

namespace {

sp::Parameters unpack(const Eigen::VectorXd& theta) {
  sp::Parameters par;
  par.log_r = theta[0];
  par.log_K = theta[1];
  par.log_q = theta[2];
  par.log_sigma_index = theta[3];
  par.logit_B0_frac = theta[4];
  return par;
}

Eigen::VectorXd pack(const sp::Parameters& par) {
  Eigen::VectorXd theta(5);
  theta[0] = par.log_r;
  theta[1] = par.log_K;
  theta[2] = par.log_q;
  theta[3] = par.log_sigma_index;
  theta[4] = par.logit_B0_frac;
  return theta;
}

double objective_with_penalty(const sp::Data& data, const Eigen::VectorXd& theta) {
  try {
    const sp::Parameters par = unpack(theta);
    const sp::Derived d = sp::evaluate_derived(data, par);

    double penalty = 0.0;

    if (d.r < 0.02 || d.r > 2.0) {
      penalty += 1e3 * std::pow(std::log(d.r / 0.35), 2.0);
    }

    if (d.K < 200.0 || d.K > 10000.0) {
      penalty += 1e3 * std::pow(std::log(d.K / 1250.0), 2.0);
    }

    if (d.sigma_index < 0.03 || d.sigma_index > 2.0) {
      penalty += 1e3 * std::pow(std::log(d.sigma_index / 0.20), 2.0);
    }

    if (d.depletion_terminal < 0.02) {
      penalty += 1e5 * std::pow(0.02 - d.depletion_terminal, 2.0);
    }

    return sp::negative_log_likelihood(data, par) + penalty;
  } catch (...) {
    return std::numeric_limits<double>::infinity();
  }
}

Eigen::VectorXd finite_difference_gradient(const sp::Data& data,
                                           const Eigen::VectorXd& theta) {
  Eigen::VectorXd grad(theta.size());

  for (int i = 0; i < theta.size(); ++i) {
    const double step = 1e-5 * (1.0 + std::abs(theta[i]));

    Eigen::VectorXd plus = theta;
    Eigen::VectorXd minus = theta;
    plus[i] += step;
    minus[i] -= step;

    const double f_plus = objective_with_penalty(data, plus);
    const double f_minus = objective_with_penalty(data, minus);

    grad[i] = (f_plus - f_minus) / (2.0 * step);
  }

  return grad;
}

class Objective {
 public:
  explicit Objective(const sp::Data& data) : data_(data) {}

  double operator()(const Eigen::VectorXd& theta, Eigen::VectorXd& grad) {
    const double f = objective_with_penalty(data_, theta);
    grad = finite_difference_gradient(data_, theta);
    return f;
  }

 private:
  const sp::Data& data_;
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

  sp::Parameters initial;
  initial.log_r = std::log(0.28);
  initial.log_K = std::log(1800.0);
  initial.log_q = std::log(0.0010);
  initial.log_sigma_index = std::log(0.25);
  initial.logit_B0_frac = std::log(0.80 / 0.20);

  Eigen::VectorXd theta = pack(initial);

  std::cout << std::fixed << std::setprecision(6);
  std::cout << "Fit surplus production example with LBFGS++\n";
  std::cout << "==========================================\n\n";
  std::cout << "initial objective = "
            << objective_with_penalty(data, theta) << "\n\n";

  LBFGSpp::LBFGSParam<double> param;
  param.epsilon = 1e-6;
  param.max_iterations = 500;
  param.max_linesearch = 50;
  param.m = 8;

  LBFGSpp::LBFGSSolver<double> solver(param);
  Objective objective(data);

  double final_objective = 0.0;
  int iterations = 0;

  try {
    iterations = solver.minimize(objective, theta, final_objective);
  } catch (const std::exception& e) {
    std::cerr << "LBFGS++ failed: " << e.what() << "\n";
    return 1;
  }

  const sp::Parameters estimated = unpack(theta);
  const Eigen::VectorXd final_grad = finite_difference_gradient(data, theta);

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

cat > run_fit_surplus_production_lbfgs_example.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/examples

set -x
"${CXX}" ${CXXFLAGS} \
  -Iexternal/Eigen \
  -Iexternal/LBFGSpp/include \
  -Iexamples/surplus_production \
  examples/surplus_production/fit_surplus_production_lbfgs.cpp \
  -o build/examples/fit_surplus_production_lbfgs

./build/examples/fit_surplus_production_lbfgs
EOF

chmod +x run_fit_surplus_production_lbfgs_example.sh

cat <<'EOF'

Installed missing LBFGS++ surplus production example.

Run:
  ./run_fit_surplus_production_lbfgs_example.sh

EOF
