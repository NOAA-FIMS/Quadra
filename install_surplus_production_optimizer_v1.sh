#!/usr/bin/env bash
set -euo pipefail

# install_surplus_production_optimizer_v1.sh
#
# Adds a simple self-contained optimizer to the surplus production example.
#
# This is intentionally not Quadra-AD yet. It gives us a fisheries-facing
# parameter-estimation example first:
#
#   estimate log_r, log_K, log_q, log_sigma_index, logit_B0_frac
#
# using finite-difference gradients + backtracking gradient descent.
#
# Next step after this works:
#   replace finite differences with Quadra gradients.

mkdir -p examples/surplus_production .quadra_patch_backups

header="examples/surplus_production/surplus_production.hpp"
if [[ ! -f "$header" ]]; then
  echo "ERROR: missing $header"
  echo "Run install_surplus_production_example_and_cleanup_v1.sh first."
  exit 1
fi

cp "$header" ".quadra_patch_backups/surplus_production.hpp.optimizer_v1.$(date +%Y%m%d_%H%M%S).bak"

cat > examples/surplus_production/fit_surplus_production.cpp <<'EOF'
#include "surplus_production.hpp"

#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

namespace sp = quadra_examples::surplus_production;

namespace {

using Theta = std::array<double, 5>;

sp::Parameters unpack(const Theta& theta) {
  sp::Parameters par;
  par.log_r = theta[0];
  par.log_K = theta[1];
  par.log_q = theta[2];
  par.log_sigma_index = theta[3];
  par.logit_B0_frac = theta[4];
  return par;
}

Theta pack(const sp::Parameters& par) {
  return Theta{
      par.log_r,
      par.log_K,
      par.log_q,
      par.log_sigma_index,
      par.logit_B0_frac};
}

double objective(const sp::Data& data, const Theta& theta) {
  try {
    const sp::Parameters par = unpack(theta);
    const sp::Derived d = sp::evaluate_derived(data, par);

    // Soft penalties keep the toy optimizer away from biologically silly or
    // numerically unstable regions.
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

Theta finite_difference_gradient(const sp::Data& data, const Theta& theta) {
  Theta grad{};

  for (std::size_t i = 0; i < theta.size(); ++i) {
    const double step = 1e-5 * (1.0 + std::abs(theta[i]));

    Theta plus = theta;
    Theta minus = theta;
    plus[i] += step;
    minus[i] -= step;

    const double f_plus = objective(data, plus);
    const double f_minus = objective(data, minus);

    grad[i] = (f_plus - f_minus) / (2.0 * step);
  }

  return grad;
}

double norm2(const Theta& x) {
  double out = 0.0;
  for (const double v : x) {
    out += v * v;
  }
  return std::sqrt(out);
}

struct FitResult {
  Theta theta;
  double objective = 0.0;
  double grad_norm = 0.0;
  int iterations = 0;
  bool converged = false;
};

FitResult fit(const sp::Data& data, const sp::Parameters& initial) {
  Theta theta = pack(initial);

  double f = objective(data, theta);
  double step_scale = 0.05;

  FitResult result;
  result.theta = theta;
  result.objective = f;

  for (int iter = 0; iter < 500; ++iter) {
    const Theta grad = finite_difference_gradient(data, theta);
    const double gnorm = norm2(grad);

    result.iterations = iter;
    result.grad_norm = gnorm;
    result.objective = f;
    result.theta = theta;

    if (iter % 25 == 0) {
      std::cout << "iter " << std::setw(4) << iter
                << " objective " << std::setw(12) << std::setprecision(6) << f
                << " grad_norm " << std::setw(12) << gnorm
                << " step " << step_scale
                << "\n";
    }

    if (gnorm < 1e-5) {
      result.converged = true;
      return result;
    }

    bool accepted = false;
    double local_step = step_scale;

    for (int ls = 0; ls < 30; ++ls) {
      Theta candidate = theta;
      for (std::size_t i = 0; i < theta.size(); ++i) {
        candidate[i] -= local_step * grad[i] / (1.0 + gnorm);
      }

      const double f_candidate = objective(data, candidate);

      if (std::isfinite(f_candidate) && f_candidate < f) {
        theta = candidate;
        f = f_candidate;
        step_scale = std::min(local_step * 1.2, 0.50);
        accepted = true;
        break;
      }

      local_step *= 0.5;
    }

    if (!accepted) {
      step_scale *= 0.25;
      if (step_scale < 1e-10) {
        result.converged = false;
        return result;
      }
    }
  }

  result.theta = theta;
  result.objective = f;
  result.grad_norm = norm2(finite_difference_gradient(data, theta));
  result.iterations = 500;
  return result;
}

void print_parameter_comparison(const sp::Parameters& initial,
                                const sp::Parameters& estimated) {
  const sp::Data data = sp::make_demo_data();
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

  std::cout << std::fixed << std::setprecision(6);

  std::cout << "Fit surplus production example\n";
  std::cout << "==============================\n\n";

  const double initial_objective = sp::negative_log_likelihood(data, initial);
  std::cout << "initial objective = " << initial_objective << "\n\n";

  const FitResult result = fit(data, initial);
  const sp::Parameters estimated = unpack(result.theta);

  std::cout << "\nFit summary\n";
  std::cout << "  converged  = " << (result.converged ? "yes" : "no") << "\n";
  std::cout << "  iterations = " << result.iterations << "\n";
  std::cout << "  objective  = " << result.objective << "\n";
  std::cout << "  grad_norm  = " << result.grad_norm << "\n";

  print_parameter_comparison(initial, estimated);

  std::cout << "\nEstimated model report\n";
  std::cout << "======================\n";
  sp::print_report(data, estimated);

  return 0;
}
EOF

cat > run_fit_surplus_production_example.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/examples

set -x
"${CXX}" ${CXXFLAGS} -Iexamples/surplus_production   examples/surplus_production/fit_surplus_production.cpp   -o build/examples/fit_surplus_production

./build/examples/fit_surplus_production
EOF

chmod +x run_fit_surplus_production_example.sh

# Append README section if not present.
if ! grep -q "Fit the demo data" examples/surplus_production/README.md; then
cat >> examples/surplus_production/README.md <<'EOF'

## Fit the demo data

A simple finite-difference optimizer example is available:

```bash
./run_fit_surplus_production_example.sh
```

This estimates:

```text
log_r
log_K
log_q
log_sigma_index
logit_B0_frac
```

The optimizer is intentionally simple. It is meant to establish the fisheries-facing
parameter-estimation workflow before replacing finite-difference gradients with
Quadra automatic differentiation.
EOF
fi

cat <<'EOF'

Installed surplus production fitting example.

Added:
  examples/surplus_production/fit_surplus_production.cpp
  run_fit_surplus_production_example.sh

Run:
  ./run_fit_surplus_production_example.sh

EOF
