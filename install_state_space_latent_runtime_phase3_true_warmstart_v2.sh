#!/usr/bin/env bash
set -euo pipefail

target="examples/state_space_surplus_production/laplace_state_space_surplus_latent_runtime.cpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  echo "Run install_state_space_latent_runtime_phase2_v2.sh first."
  exit 1
fi

mkdir -p .quadra_patch_backups
cp "$target" ".quadra_patch_backups/laplace_state_space_surplus_latent_runtime.cpp.phase3_warmstart_v2.$(date +%Y%m%d_%H%M%S).bak"

python3 - "$target" <<'PY'
import re
import sys
from pathlib import Path

p = Path(sys.argv[1])
s = p.read_text()

if 'optimize_x_from(const ss::Data& data' in s:
    print('Phase 3 warm-start patch already appears to be installed.')
    raise SystemExit(0)

old_sig = 'Eigen::VectorXd optimize_x(const ss::Data& data, const ss::Parameters& par) {'
new_sig = 'Eigen::VectorXd optimize_x_from(const ss::Data& data, const ss::Parameters& par, const Eigen::VectorXd& initial_x) {'

sig_pos = s.find(old_sig)
if sig_pos < 0:
    raise SystemExit('Could not find optimize_x signature')

s = s[:sig_pos] + new_sig + s[sig_pos + len(old_sig):]

func_start = sig_pos
next_markers = [
    '\nEigen::SparseMatrix',
    '\nstruct EvalResult',
    '\nclass LatentTridiagonalLaplaceRuntime',
]
func_end_candidates = [s.find(m, func_start + len(new_sig)) for m in next_markers]
func_end_candidates = [x for x in func_end_candidates if x >= 0]
func_end = min(func_end_candidates) if func_end_candidates else len(s)
func = s[func_start:func_end]

pattern = re.compile(r'  Eigen::VectorXd x\s*=\s*[^;]+;')
m = pattern.search(func)
if not m:
    raise SystemExit('Could not find local Eigen::VectorXd x initialization inside optimize_x')

func = func[:m.start()] + '  Eigen::VectorXd x = initial_x;' + func[m.end():]
s = s[:func_start] + func + s[func_end:]

wrapper_marker = 'class LatentTridiagonalLaplaceRuntime {'
if wrapper_marker not in s:
    raise SystemExit('Could not find runtime class marker')

wrapper = (
    '\nEigen::VectorXd optimize_x(const ss::Data& data, const ss::Parameters& par) {\n'
    '  return optimize_x_from(\n'
    '      data,\n'
    '      par,\n'
    '      Eigen::VectorXd::Zero(\n'
    '          static_cast<int>(data.index_observed.size() - 1)));\n'
    '}\n\n'
)
s = s.replace(wrapper_marker, wrapper + wrapper_marker, 1)

old_call = '    xhat_ = optimize_x(data_, par_);'
new_call = '    xhat_ = optimize_x_from(data_, par_, xhat_);'

idx1 = s.find(old_call)
if idx1 < 0:
    raise SystemExit('Could not find cold optimize_x runtime call')

idx2 = s.find(old_call, idx1 + len(old_call))
if idx2 < 0:
    raise SystemExit('Could not find warm optimize_x runtime call')

s = s[:idx2] + new_call + s[idx2 + len(old_call):]

s = s.replace(
    '    // Phase 2 persists xhat/backend first. Phase 3 should expose\n'
    '    // optimize_x_from(start) so this becomes a true warm start.\n'
    '    xhat_ = optimize_x_from(data_, par_, xhat_);',
    '    // True warm start from previously cached xhat_.\n'
    '    xhat_ = optimize_x_from(data_, par_, xhat_);'
)

p.write_text(s)
print('Installed true warm-start patch.')
PY

cat <<'EOF'

Installed Phase 3 true warm-start runtime patch.

Run:
  ./run_state_space_surplus_latent_runtime_phase2.sh 20

Expected:
  backend = tridiagonal
  objective = -10.642176

EOF
