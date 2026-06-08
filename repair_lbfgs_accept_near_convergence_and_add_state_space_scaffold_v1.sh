#!/usr/bin/env bash
set -euo pipefail

mkdir -p .quadra_patch_backups examples/state_space_surplus_production

target="examples/surplus_production/fit_surplus_production_lbfgs.cpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  echo "Run add_missing_surplus_production_lbfgs_example_v1.sh first."
  exit 1
fi

cp "$target" ".quadra_patch_backups/fit_surplus_production_lbfgs.cpp.accept_near_convergence.$(date +%Y%m%d_%H%M%S).bak"

python3 - <<'PYEOF'
from pathlib import Path

p = Path("examples/surplus_production/fit_surplus_production_lbfgs.cpp")
s = p.read_text()

old = '''  double final_objective = 0.0;
  int iterations = 0;

  try {
    iterations = solver.minimize(objective, x, final_objective);
  } catch (const std::exception& e) {
    std::cerr << "LBFGS++ failed: " << e.what() << "\\n";
    std::cerr << "Current objective = "
              << objective_value(data, center, x) << "\\n";
    std::cerr << "Current gradient norm = "
              << finite_difference_gradient(data, center, x).norm() << "\\n";
    return 1;
  }

  const sp::Parameters estimated = unpack_scaled(x, center);
  const Eigen::VectorXd final_grad = finite_difference_gradient(data, center, x);

  std::cout << "Fit summary\\n";
  std::cout << "  iterations = " << iterations << "\\n";
  std::cout << "  objective  = " << final_objective << "\\n";
  std::cout << "  grad_norm  = " << final_grad.norm() << "\\n";
'''

new = '''  double final_objective = 0.0;
  int iterations = 0;
  bool converged = false;
  bool accepted_line_search_failure = false;

  try {
    iterations = solver.minimize(objective, x, final_objective);
    converged = true;
  } catch (const std::exception& e) {
    final_objective = objective_value(data, center, x);
    const Eigen::VectorXd current_grad =
        finite_difference_gradient(data, center, x);
    const double current_grad_norm = current_grad.norm();

    if (std::isfinite(final_objective) && current_grad_norm < 1e-2) {
      accepted_line_search_failure = true;
      converged = true;

      std::cout << "LBFGS++ terminated during line search after reaching "
                << "an acceptable finite-difference optimum.\\n";
      std::cout << "  message    = " << e.what() << "\\n";
      std::cout << "  objective  = " << final_objective << "\\n";
      std::cout << "  grad_norm  = " << current_grad_norm << "\\n\\n";
    } else {
      std::cerr << "LBFGS++ failed: " << e.what() << "\\n";
      std::cerr << "Current objective = " << final_objective << "\\n";
      std::cerr << "Current gradient norm = " << current_grad_norm << "\\n";
      return 1;
    }
  }

  const sp::Parameters estimated = unpack_scaled(x, center);
  const Eigen::VectorXd final_grad = finite_difference_gradient(data, center, x);

  std::cout << "Fit summary\\n";
  std::cout << "  converged  = " << (converged ? "yes" : "no") << "\\n";
  std::cout << "  accepted line-search termination = "
            << (accepted_line_search_failure ? "yes" : "no") << "\\n";
  std::cout << "  iterations = " << iterations << "\\n";
  std::cout << "  objective  = " << final_objective << "\\n";
  std::cout << "  grad_norm  = " << final_grad.norm() << "\\n";
'''

if old not in s:
    raise SystemExit("Could not find LBFGS++ try/catch block to patch")

s = s.replace(old, new, 1)
p.write_text(s)
PYEOF

cat > examples/state_space_surplus_production/README.md <<'EOF'
# State-space surplus production scaffold

This is the next fisheries-facing Quadra example after the deterministic
Schaefer surplus production model.

The intended model is:

```text
B[t+1] = B[t] + r B[t] (1 - B[t] / K) - C[t] + process error
I[t]   = q B[t] exp(observation error)
```

A numerically safer implementation will usually work on log biomass:

```text
log_B[t+1] = log(predicted_B[t+1]) + eta[t]
eta[t] ~ Normal(0, sigma_process)
```

with index observations:

```text
log(I[t]) ~ Normal(log(q) + log_B[t], sigma_index)
```

## Fixed effects

A first version should estimate:

```text
log_r
log_K
log_q
log_sigma_process
log_sigma_index
logit_B0_frac
```

## Random effects

The random effects are annual biomass/process deviations:

```text
eta[0], eta[1], ..., eta[n - 2]
```

or equivalently latent annual log biomass states after the initial condition.

## Laplace target

The objective should be separable into:

```text
joint_nll(theta, u)
```

where:

```text
theta = fixed effects
u     = latent process deviations / biomass states
```

Then Quadra's Laplace machinery can evaluate:

```text
marginal_nll(theta)
=
joint_nll(theta, u_hat)
+
0.5 log |H_uu(theta, u_hat)|
-
n_u/2 log(2 pi)
```
EOF

cat > examples/state_space_surplus_production/run_state_space_surplus_production.cpp <<'EOF'
#include <iostream>

int main() {
  std::cout << "State-space surplus production scaffold\n";
  std::cout << "Next step: implement joint_objective(theta, u) with latent biomass/process effects.\n";
  return 0;
}
EOF

cat > run_state_space_surplus_production_scaffold.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/examples

set -x
"${CXX}" ${CXXFLAGS}   examples/state_space_surplus_production/run_state_space_surplus_production.cpp   -o build/examples/run_state_space_surplus_production

./build/examples/run_state_space_surplus_production
EOF

chmod +x run_state_space_surplus_production_scaffold.sh

cat <<'EOF'

Patched LBFGS++ near-convergence handling and added state-space surplus production scaffold.

Run:
  ./run_fit_surplus_production_lbfgs_example.sh
  ./run_state_space_surplus_production_scaffold.sh

EOF
