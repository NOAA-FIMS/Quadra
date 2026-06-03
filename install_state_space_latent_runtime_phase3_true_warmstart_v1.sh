#!/usr/bin/env bash
set -euo pipefail

target="examples/state_space_surplus_production/laplace_state_space_surplus_latent_runtime.cpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  echo "Run install_state_space_latent_runtime_phase2_v2.sh first."
  exit 1
fi

mkdir -p .quadra_patch_backups
cp "$target" ".quadra_patch_backups/laplace_state_space_surplus_latent_runtime.cpp.phase3_warmstart.$(date +%Y%m%d_%H%M%S).bak"

python3 - "$target" <<'PY'
import sys
from pathlib import Path

p = Path(sys.argv[1])
s = p.read_text()

if "optimize_x_from(const ss::Data& data" in s:
    print("Phase 3 warm-start patch already appears to be installed.")
    raise SystemExit(0)

old_sig = "Eigen::VectorXd optimize_x(const ss::Data& data, const ss::Parameters& par) {"
new_sig = "Eigen::VectorXd optimize_x_from(const ss::Data& data, const ss::Parameters& par, const Eigen::VectorXd& initial_x) {"

if old_sig not in s:
    raise SystemExit("Could not find optimize_x signature")

s = s.replace(old_sig, new_sig, 1)

old_init_candidates = [
    "  Eigen::VectorXd x = Eigen::VectorXd::Zero(static_cast<int>(data.index_observed.size() - 1));",
    "  Eigen::VectorXd x = Eigen::VectorXd::Zero(data.index_observed.size() - 1);",
]

for old_init in old_init_candidates:
    if old_init in s:
        s = s.replace(old_init, "  Eigen::VectorXd x = initial_x;", 1)
        break
else:
    raise SystemExit("Could not find optimize_x zero initialization")

wrapper_marker = "class LatentTridiagonalLaplaceRuntime {"
if wrapper_marker not in s:
    raise SystemExit("Could not find runtime class marker")

wrapper = '''
Eigen::VectorXd optimize_x(const ss::Data& data, const ss::Parameters& par) {
  return optimize_x_from(
      data,
      par,
      Eigen::VectorXd::Zero(
          static_cast<int>(data.index_observed.size() - 1)));
}

'''

s = s.replace(wrapper_marker, wrapper + wrapper_marker, 1)

old_warm = "    xhat_ = optimize_x(data_, par_);"
new_warm = "    xhat_ = optimize_x_from(data_, par_, xhat_);"

idx1 = s.find(old_warm)
if idx1 < 0:
    raise SystemExit("Could not find first optimize_x runtime call")

idx2 = s.find(old_warm, idx1 + len(old_warm))
if idx2 < 0:
    raise SystemExit("Could not find second optimize_x runtime call for warm path")

s = s[:idx2] + new_warm + s[idx2 + len(old_warm):]

s = s.replace(
    "// Phase 2 persists xhat/backend first. Phase 3 should expose\n"
    "    // optimize_x_from(start) so this becomes a true warm start.\n"
    "    xhat_ = optimize_x_from(data_, par_, xhat_);",
    "// True warm start from previously cached xhat_.\n"
    "    xhat_ = optimize_x_from(data_, par_, xhat_);"
)

p.write_text(s)
PY

cat <<'EOF'

Installed Phase 3 true warm-start runtime patch.

Run:
  ./run_state_space_surplus_latent_runtime_phase2.sh 20

Expected:
  backend = tridiagonal
  objective = -10.642176

Watch:
  warm_avg_ms should drop relative to Phase 2 if the inner optimizer benefits from xhat_.

EOF
