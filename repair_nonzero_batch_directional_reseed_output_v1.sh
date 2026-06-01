#!/usr/bin/env bash
set -euo pipefail

mkdir -p .quadra_patch_backups

target="tests/test_had_quadra_nonzero_batch_directional.cpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/test_had_quadra_nonzero_batch_directional.cpp.reseed_output.$(date +%Y%m%d_%H%M%S).bak"

python3 - <<'PYEOF'
from pathlib import Path

p = Path("tests/test_had_quadra_nonzero_batch_directional.cpp")
s = p.read_text()

old = "    had::PropagateAdjointDirectionalBatch();"
new = """    // PropagateAdjoint() consumes/clears reverse adjoints while building the
    // base Hessian, so the directional reverse pass needs the output adjoint
    // seed restored.
    had::g_ADGraph->vertices[f.varId].w = 1.0;
    had::PropagateAdjointDirectionalBatch();"""

if old not in s:
    raise SystemExit("Could not find PropagateAdjointDirectionalBatch call")

s = s.replace(old, new, 1)
p.write_text(s)
PYEOF

cat <<'EOF'

Patched nonzero diagnostic to reseed output adjoint before batch reverse.

Run:
  ./run_had_quadra_nonzero_batch_directional_test.sh

EOF
