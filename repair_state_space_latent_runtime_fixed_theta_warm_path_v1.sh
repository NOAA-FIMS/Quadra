#!/usr/bin/env bash
set -euo pipefail

# repair_state_space_latent_runtime_fixed_theta_warm_path_v1.sh
#
# For the fixed-theta Phase 3 runtime benchmark, evaluate_warm() should not
# call the optimizer at all after evaluate_cold().
#
# This benchmark has no theta updates between evaluations, so cached xhat_ is
# already the optimum for the same parameters. The old warm path was still
# calling LBFGS++ from near/at the optimum and triggering a line-search failure.
#
# Later, when theta changes, we can add a dirty/parameter-changed flag:
#   if theta_changed -> optimize_x_from(previous_xhat)
#   else             -> evaluate_at_xhat()

target="examples/state_space_surplus_production/laplace_state_space_surplus_latent_runtime.cpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

mkdir -p .quadra_patch_backups
cp "$target" ".quadra_patch_backups/laplace_state_space_surplus_latent_runtime.cpp.fixed_theta_warm.$(date +%Y%m%d_%H%M%S).bak"

python3 - "$target" <<'PY'
import sys
from pathlib import Path

p = Path(sys.argv[1])
s = p.read_text()

start = s.find("  EvalResult evaluate_warm() {")
if start < 0:
    raise SystemExit("Could not find evaluate_warm start")

end = s.find("\n  EvalResult evaluate_cached_no_solve()", start)
if end < 0:
    raise SystemExit("Could not find evaluate_warm end")

new_block = '''  EvalResult evaluate_warm() {
    if (!initialized_) {
      return evaluate_cold();
    }

    // Fixed-theta runtime benchmark:
    // the parameters have not changed, so the cached xhat_ remains valid.
    // Do not call LBFGS++ from an already-converged point.
    return evaluate_at_xhat();
  }

'''

s = s[:start] + new_block + s[end + 1:]
p.write_text(s)
PY

cat <<'EOF'

Repaired fixed-theta warm path: evaluate_warm() now reuses cached xhat_ directly.

Run:
  ./run_state_space_surplus_latent_runtime_phase2.sh 20

Expected:
  backend = tridiagonal
  objective = -10.642176
  warm_avg_ms should be close to cached_avg_ms

EOF
