#!/usr/bin/env bash
set -euo pipefail

# repair_state_space_latent_runtime_cold_initializer_v1.sh
#
# Restores the model-aware cold initializer after Phase 3.
#
# Correct behavior:
#   optimize_x_from(data, par, initial_x)  uses supplied initial_x
#   optimize_x(data, par)                 starts from deterministic_initial_x(data, par)
#   evaluate_warm()                       reuses cached xhat_, or in fixed-theta mode skips solve

target="examples/state_space_surplus_production/laplace_state_space_surplus_latent_runtime.cpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

mkdir -p .quadra_patch_backups
cp "$target" ".quadra_patch_backups/laplace_state_space_surplus_latent_runtime.cpp.cold_initializer.$(date +%Y%m%d_%H%M%S).bak"

python3 - "$target" <<'PY'
import sys
from pathlib import Path

p = Path(sys.argv[1])
s = p.read_text()

old = '''Eigen::VectorXd optimize_x(const ss::Data& data, const ss::Parameters& par) {
  return optimize_x_from(
      data,
      par,
      Eigen::VectorXd::Zero(
          static_cast<int>(data.index_observed.size() - 1)));
}'''

new = '''Eigen::VectorXd optimize_x(const ss::Data& data, const ss::Parameters& par) {
  return optimize_x_from(
      data,
      par,
      deterministic_initial_x(data, par));
}'''

if old not in s:
    raise SystemExit("Could not find zero-start optimize_x wrapper")

s = s.replace(old, new, 1)
p.write_text(s)
PY

cat <<'EOF'

Restored deterministic cold initializer.

Run:
  ./run_state_space_surplus_latent_runtime_phase2.sh 20

Expected:
  backend = tridiagonal
  objective = -10.642176
  cold solve should no longer fail.

EOF
