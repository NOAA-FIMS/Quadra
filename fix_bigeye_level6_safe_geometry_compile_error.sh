#!/usr/bin/env bash
set -euo pipefail

L6="examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity"
SAFE="$L6/diagnostics/bigeye_safe_fixed_effect_diagnostics.hpp"
DRIVER="$L6/quadra/bigeye_level6_purse_seine_age_selectivity.cpp"

if [[ ! -f "$SAFE" || ! -f "$DRIVER" ]]; then
  echo "ERROR: expected Level 6 files not found. Run from repo root."
  exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
cp "$SAFE" "${SAFE}.before_safe_geometry_compile_fix.${STAMP}"
cp "$DRIVER" "${DRIVER}.before_safe_geometry_compile_fix.${STAMP}"

python3 - <<'PY'
from pathlib import Path
import re

safe = Path("examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity/diagnostics/bigeye_safe_fixed_effect_diagnostics.hpp")
s = safe.read_text()

start = s.find("// Lightweight safe fixed-effect geometry summary.")
if start != -1:
    s = s[:start].rstrip() + "\n"

safe.write_text(s)

driver = Path("examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity/quadra/bigeye_level6_purse_seine_age_selectivity.cpp")
d = driver.read_text()

d = re.sub(
    r'\n\s*pifsc_bigeye_tuna::write_safe_fixed_effect_geometry_summary\(\n'
    r'\s*"examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity/outputs/"\n'
    r'\s*"bigeye_level6_safe_fixed_effect_geometry_summary\.txt",\n'
    r'\s*"examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity/outputs/"\n'
    r'\s*"bigeye_level6_safe_fixed_effect_geometry_summary\.csv",\n'
    r'\s*params, fit\);\n',
    "\n",
    d,
)

d = re.sub(
    r'\n\s*std::cout << "wrote:\s+"\n'
    r'\s*<< "examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity/outputs/"\n'
    r'\s*"bigeye_level6_safe_fixed_effect_geometry_summary\.txt\\n";\n'
    r'\s*std::cout << "wrote:\s+"\n'
    r'\s*<< "examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity/outputs/"\n'
    r'\s*"bigeye_level6_safe_fixed_effect_geometry_summary\.csv\\n";\n',
    "\n",
    d,
)

d = d.replace('\\\\n";', '\\n";')
driver.write_text(d)
PY

cat > inspect_bigeye_level6_compile_fix.sh <<'EOF_INSPECT'
#!/usr/bin/env bash
set -euo pipefail

L6="examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity"
SAFE="$L6/diagnostics/bigeye_safe_fixed_effect_diagnostics.hpp"
DRIVER="$L6/quadra/bigeye_level6_purse_seine_age_selectivity.cpp"

echo "== Removed broken safe geometry references? =="
grep -n "SafeFixedEffectGeometry\\|safe_fixed_effect_geometry\\|fixed_hessian" "$SAFE" "$DRIVER" || true

echo
echo "== Safe wiggle references still present? =="
grep -n "safe_fixed_effect_wiggle\\|write_safe_fixed_effect_wiggle" "$SAFE" "$DRIVER" | head -80

echo
echo "== Literal escaped newline check =="
grep -n '\\\\n' "$DRIVER" || true
EOF_INSPECT

chmod +x inspect_bigeye_level6_compile_fix.sh

cat > run_bigeye_level6_compile_fix_check.sh <<'EOF_RUN'
#!/usr/bin/env bash
set -euo pipefail

./inspect_bigeye_level6_compile_fix.sh

echo
echo "== O3 build Bigeye Level 6 after safe-geometry compile fix =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity/quadra/bigeye_level6_purse_seine_age_selectivity.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level6_compile_fix_check

echo
echo "== Run Bigeye Level 6 after safe-geometry compile fix =="
./build/examples/pifsc_bigeye_level6_compile_fix_check

echo
echo "== Safe wiggle diagnostics preview =="
sed -n '1,140p' \
  examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity/outputs/bigeye_level6_safe_fixed_effect_wiggle_diagnostics.txt

echo
echo "== Commit-ready status =="
git status --short examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity examples/NMFS/pifsc_bigeye_tuna/workflow/scientific_reasoning_log.md
EOF_RUN

chmod +x run_bigeye_level6_compile_fix_check.sh

echo "Patched Level 6 compile failure by removing broken safe geometry block."
echo
echo "Run:"
echo "  ./inspect_bigeye_level6_compile_fix.sh"
echo "  ./run_bigeye_level6_compile_fix_check.sh"
