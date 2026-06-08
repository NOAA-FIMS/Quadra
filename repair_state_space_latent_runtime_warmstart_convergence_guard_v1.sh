#!/usr/bin/env bash
set -euo pipefail

# repair_state_space_latent_runtime_warmstart_convergence_guard_v1.sh
#
# Fixes Phase 3 true warm-start failure when xhat_ is already at/near the optimum.
#
# Problem:
#   optimize_x_from(data, par, xhat_) can throw an LBFGS++ line-search error
#   when started essentially at the optimum.
#
# Fix:
#   In evaluate_warm(), check fd_grad_x(data_, par_, xhat_).norm().
#   If already small, skip the optimizer and directly evaluate at cached xhat_.
#   Also catch optimizer exceptions and accept cached xhat_ if gradient is small.

target="examples/state_space_surplus_production/laplace_state_space_surplus_latent_runtime.cpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

mkdir -p .quadra_patch_backups
cp "$target" ".quadra_patch_backups/laplace_state_space_surplus_latent_runtime.cpp.warm_guard.$(date +%Y%m%d_%H%M%S).bak"

python3 - "$target" <<'PY'
import sys
from pathlib import Path

p = Path(sys.argv[1])
s = p.read_text()

old = '''  EvalResult evaluate_warm() {
    if (!initialized_) {
      return evaluate_cold();
    }

    // True warm start from previously cached xhat_.
    xhat_ = optimize_x_from(data_, par_, xhat_);
    initialized_ = true;
    return evaluate_at_xhat();
  }'''

new = '''  EvalResult evaluate_warm() {
    if (!initialized_) {
      return evaluate_cold();
    }

    const double cached_grad_norm =
        fd_grad_x(data_, par_, xhat_).norm();

    if (cached_grad_norm < 1e-5) {
      return evaluate_at_xhat();
    }

    try {
      xhat_ = optimize_x_from(data_, par_, xhat_);
    } catch (const std::exception&) {
      const double fallback_grad_norm =
          fd_grad_x(data_, par_, xhat_).norm();

      if (fallback_grad_norm < 1e-5) {
        return evaluate_at_xhat();
      }

      throw;
    }

    initialized_ = true;
    return evaluate_at_xhat();
  }'''

if old not in s:
    raise SystemExit("Could not find exact evaluate_warm block")

s = s.replace(old, new, 1)
p.write_text(s)
PY

cat <<'EOF'

Installed warm-start convergence guard.

Run:
  ./run_state_space_surplus_latent_runtime_phase2.sh 20

Expected:
  backend = tridiagonal
  objective = -10.642176
  warm_avg_ms should now be close to cached_avg_ms

EOF
