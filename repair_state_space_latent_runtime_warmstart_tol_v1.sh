#!/usr/bin/env bash
set -euo pipefail

# repair_state_space_latent_runtime_warmstart_tol_v1.sh
#
# The cached xhat has grad_norm around 9.8e-5, but the previous guard used 1e-5.
# That still calls LBFGS++ from an effectively converged point and can trigger
# a line-search failure. Use 1e-3, consistent with prior accepted solve logic.

target="examples/state_space_surplus_production/laplace_state_space_surplus_latent_runtime.cpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

mkdir -p .quadra_patch_backups
cp "$target" ".quadra_patch_backups/laplace_state_space_surplus_latent_runtime.cpp.warm_tol.$(date +%Y%m%d_%H%M%S).bak"

python3 - "$target" <<'PY'
import sys
from pathlib import Path

p = Path(sys.argv[1])
s = p.read_text()

count = s.count("< 1e-5")
if count == 0:
    raise SystemExit("Did not find any '< 1e-5' warm-start guards to relax")

s = s.replace("< 1e-5", "< 1e-3")
p.write_text(s)

print(f"Relaxed {count} warm-start convergence guard(s): 1e-5 -> 1e-3")
PY

cat <<'EOF'

Relaxed warm-start convergence tolerance.

Run:
  ./run_state_space_surplus_latent_runtime_phase2.sh 20

Expected:
  backend = tridiagonal
  objective = -10.642176
  warm_avg_ms should be near cached_avg_ms because cached xhat is accepted.

EOF
