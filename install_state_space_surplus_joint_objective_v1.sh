#!/usr/bin/env bash
set -euo pipefail

# install_state_space_surplus_joint_objective_v1.sh
#
# Adds the first concrete state-space surplus production example:
#   - joint_objective(theta, u)
#   - fixed effects theta
#   - random effects/process deviations u
#   - log-biomass reconstruction
#   - process likelihood
#   - index observation likelihood
#   - fisheries reference point report
#
# This is not Laplace yet. It is the joint objective foundation:
#
#   joint_nll(theta, u)
#
# Next step:
#   solve for u_hat(theta), then add Laplace correction.

mkdir -p examples/state_space_surplus_production .quadra_patch_backups

if [[ ! -f examples/surplus_production/surplus_production.hpp ]]; then
  echo "ERROR: missing examples/surplus_production/surplus_production.hpp"
  echo "Run the surplus production example installer first."
  exit 1
fi

cat > examples/state_space_surplus_production/state_space_surplus_production.hpp <<'EOF'
#ifndef QUADRA_EXAMPLES_STATE_SPACE_SURPLUS_PRODUCTION_HPP
#define QUADRA_EXAMPLES_STATE_SPACE_SURPLUS_PRODUCTION_HPP

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../surplus_production/surplus_production.hpp"

namespace quadra_examples {
namespace state_space_surplus_production {

struct Data {
  std::vector<double> catch_observed;
  std::vector<double> index_observed;
};

struct Parameters {
  double log_r = std::log(0.45);
  double log_K = std::log(900.0);
  double log_q = std::log(0.0020);
  double log_sigma_process = std::log(0.12);
  double log_sigma_index = std::log(0.12);
  double logit_B0_frac = std::log(0.85 / 0.15);
};

struct Derived {
  double r = 0.0;
  double K = 0.0;
  double q = 0.0;
  double sigma_process = 0.0;
  double sigma_index = 0.0;
  double B0_frac = 0.0;
  double MSY = 0.0;
  double B_MSY = 0.0;
  double F_MSY = 0.0;
  double depletion_terminal = 0.0;
  std::vector<double> log_biomass;
  std::vector<double> biomass;
  std::vector<double> deterministic_next_biomass;
  std::vector<double> process_residuals;
  std::vector<double> index_predicted;
  std::vector<double> log_index_residuals;
};

inline double inv_logit(const double x) {
  return quadra_examples::surplus_production::inv_logit(x);
}

inline Data make_demo_data() {
  const auto base = quadra_examples::surplus_production::make_demo_data();
  Data data;
  data.catch_observed = base.catch_observed;
  data.index_observed = base.index_observed;
  return data;
}

inline Parameters make_demo_parameters() {
  Parameters par;
  par.log_r = std::log(0.50);
  par.log_K = std::log(700.0);
  par.log_q = std::log(0.0024);
  par.log_sigma_process = std::log(0.15);
  par.log_sigma_index = std::log(0.10);
  par.logit_B0_frac = std::log(0.90 / 0.10);
  return par;
}

inline void validate_data(const Data& data) {
  if (data.catch_observed.empty()) {
    throw std::runtime_error("catch_observed must not be empty");
  }

  if (data.index_observed.size() != data.catch_observed.size()) {
    throw std::runtime_error(
        "index_observed and catch_observed must have same length");
  }

  for (const double c : data.catch_observed) {
    if (!(c >= 0.0)) {
      throw std::runtime_error("catch_observed must be nonnegative");
    }
  }

  for (const double i : data.index_observed) {
    if (!(i > 0.0)) {
      throw std::runtime_error("index_observed must be positive");
    }
  }
}

inline double normal_nll(const double residual, const double sigma) {
  const double z = residual / sigma;
  return std::log(sigma) + 0.5 * std::log(2.0 * M_PI) + 0.5 * z * z;
}

inline Derived evaluate_derived(const Data& data,
                                const Parameters& par,
                                const std::vector<double>& u) {
  validate_data(data);

  const int n = static_cast<int>(data.catch_observed.size());
  const int expected_u = std::max(0, n - 1);

  if (static_cast<int>(u.size()) != expected_u) {
    throw std::runtime_error("u must have length n_years - 1");
  }

  Derived out;

  out.r = std::exp(par.log_r);
  out.K = std::exp(par.log_K);
  out.q = std::exp(par.log_q);
  out.sigma_process = std::exp(par.log_sigma_process);
  out.sigma_index = std::exp(par.log_sigma_index);
  out.B0_frac = inv_logit(par.logit_B0_frac);

  if (!(out.r > 0.0) || !(out.K > 0.0) || !(out.q > 0.0) ||
      !(out.sigma_process > 0.0) || !(out.sigma_index > 0.0) ||
      !(out.B0_frac > 0.0) || !(out.B0_frac < 1.0)) {
    throw std::runtime_error("invalid transformed parameter");
  }

  out.log_biomass.assign(static_cast<std::size_t>(n), 0.0);
  out.biomass.assign(static_cast<std::size_t>(n), 0.0);
  out.deterministic_next_biomass.assign(static_cast<std::size_t>(expected_u), 0.0);
  out.process_residuals.assign(static_cast<std::size_t>(expected_u), 0.0);
  out.index_predicted.assign(static_cast<std::size_t>(n), 0.0);
  out.log_index_residuals.assign(static_cast<std::size_t>(n), 0.0);

  out.biomass[0] = out.B0_frac * out.K;
  out.log_biomass[0] = std::log(out.biomass[0]);

  for (int t = 0; t < n; ++t) {
    const double B_t = std::exp(out.log_biomass[static_cast<std::size_t>(t)]);

    out.biomass[static_cast<std::size_t>(t)] = B_t;
    out.index_predicted[static_cast<std::size_t>(t)] = out.q * B_t;
    out.log_index_residuals[static_cast<std::size_t>(t)] =
        std::log(data.index_observed[static_cast<std::size_t>(t)]) -
        std::log(out.index_predicted[static_cast<std::size_t>(t)]);

    if (t < n - 1) {
      const double production = out.r * B_t * (1.0 - B_t / out.K);
      const double deterministic_next =
          std::max(B_t + production -
                       data.catch_observed[static_cast<std::size_t>(t)],
                   1e-9);

      out.deterministic_next_biomass[static_cast<std::size_t>(t)] =
          deterministic_next;

      // Random effect u[t] is a log-scale process residual:
      //
      //   log_B[t+1] = log(deterministic_next) + u[t]
      //
      out.process_residuals[static_cast<std::size_t>(t)] = u[static_cast<std::size_t>(t)];
      out.log_biomass[static_cast<std::size_t>(t + 1)] =
          std::log(deterministic_next) + u[static_cast<std::size_t>(t)];
    }
  }

  out.MSY = out.r * out.K / 4.0;
  out.B_MSY = out.K / 2.0;
  out.F_MSY = out.r / 2.0;
  out.depletion_terminal = out.biomass.back() / out.K;

  return out;
}

inline double joint_objective(const Data& data,
                              const Parameters& par,
                              const std::vector<double>& u) {
  const Derived d = evaluate_derived(data, par, u);

  double nll = 0.0;

  for (const double residual : d.process_residuals) {
    nll += normal_nll(residual, d.sigma_process);
  }

  for (const double residual : d.log_index_residuals) {
    nll += normal_nll(residual, d.sigma_index);
  }

  return nll;
}

inline std::vector<double> zero_random_effects(const Data& data) {
  const int n = static_cast<int>(data.catch_observed.size());
  return std::vector<double>(static_cast<std::size_t>(std::max(0, n - 1)), 0.0);
}

inline void print_report(const Data& data,
                         const Parameters& par,
                         const std::vector<double>& u) {
  const Derived d = evaluate_derived(data, par, u);
  const double joint_nll = joint_objective(data, par, u);

  std::cout << std::fixed << std::setprecision(6);

  std::cout << "State-space surplus production joint objective\n";
  std::cout << "----------------------------------------------\n";
  std::cout << "joint negative log likelihood = " << joint_nll << "\n\n";

  std::cout << "Parameters\n";
  std::cout << "  r              = " << d.r << "\n";
  std::cout << "  K              = " << d.K << "\n";
  std::cout << "  q              = " << d.q << "\n";
  std::cout << "  sigma_process  = " << d.sigma_process << "\n";
  std::cout << "  sigma_index    = " << d.sigma_index << "\n";
  std::cout << "  B0 / K         = " << d.B0_frac << "\n\n";

  std::cout << "Reference points\n";
  std::cout << "  MSY            = " << d.MSY << "\n";
  std::cout << "  B_MSY          = " << d.B_MSY << "\n";
  std::cout << "  F_MSY          = " << d.F_MSY << "\n";
  std::cout << "  B_terminal / K = " << d.depletion_terminal << "\n\n";

  std::cout << std::setw(6) << "year"
            << std::setw(14) << "catch"
            << std::setw(14) << "biomass"
            << std::setw(14) << "index obs"
            << std::setw(14) << "index pred"
            << std::setw(14) << "idx resid"
            << std::setw(14) << "proc resid"
            << "\n";

  for (std::size_t t = 0; t < data.catch_observed.size(); ++t) {
    std::cout << std::setw(6) << t
              << std::setw(14) << data.catch_observed[t]
              << std::setw(14) << d.biomass[t]
              << std::setw(14) << data.index_observed[t]
              << std::setw(14) << d.index_predicted[t]
              << std::setw(14) << d.log_index_residuals[t];

    if (t < d.process_residuals.size()) {
      std::cout << std::setw(14) << d.process_residuals[t];
    } else {
      std::cout << std::setw(14) << "-";
    }

    std::cout << "\n";
  }
}

}  // namespace state_space_surplus_production
}  // namespace quadra_examples

#endif
EOF

cat > examples/state_space_surplus_production/run_state_space_surplus_production.cpp <<'EOF'
#include "state_space_surplus_production.hpp"

#include <cmath>
#include <iostream>

namespace ss = quadra_examples::state_space_surplus_production;

int main() {
  const ss::Data data = ss::make_demo_data();
  const ss::Parameters par = ss::make_demo_parameters();
  const std::vector<double> u = ss::zero_random_effects(data);

  ss::print_report(data, par, u);

  return 0;
}
EOF

cat > run_state_space_surplus_production_joint_example.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/examples

set -x
"${CXX}" ${CXXFLAGS} \
  -Iexamples/state_space_surplus_production \
  -Iexamples/surplus_production \
  examples/state_space_surplus_production/run_state_space_surplus_production.cpp \
  -o build/examples/run_state_space_surplus_production_joint

./build/examples/run_state_space_surplus_production_joint
EOF

chmod +x run_state_space_surplus_production_joint_example.sh

cat > examples/state_space_surplus_production/README.md <<'EOF'
# State-space surplus production example

This example extends the deterministic Schaefer surplus production model with
log-scale process deviations in biomass dynamics.

The joint objective is:

```text
joint_nll(theta, u)
```

where:

```text
theta = fixed effects
u     = annual log-scale process deviations
```

The process model is:

```text
pred_B[t+1] = B[t] + r B[t] (1 - B[t] / K) - C[t]
log_B[t+1]  = log(pred_B[t+1]) + u[t]
u[t]        ~ Normal(0, sigma_process)
```

The observation model is:

```text
log(I[t]) ~ Normal(log(q) + log_B[t], sigma_index)
```

Run:

```bash
./run_state_space_surplus_production_joint_example.sh
```

Next steps:

1. Add finite-difference gradient checks for `joint_objective(theta, u)`.
2. Add a Newton solver for `u_hat(theta)`.
3. Add the Laplace correction.
4. Replace finite differences with Quadra exact gradients.
EOF

cat <<'EOF'

Installed state-space surplus production joint objective.

Added:
  examples/state_space_surplus_production/state_space_surplus_production.hpp
  examples/state_space_surplus_production/run_state_space_surplus_production.cpp
  run_state_space_surplus_production_joint_example.sh

Run:
  ./run_state_space_surplus_production_joint_example.sh

EOF
