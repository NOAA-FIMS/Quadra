#!/usr/bin/env bash
set -euo pipefail

DRIVER="examples/NMFS/pifsc_bigeye_tuna/level10_q_anchor/quadra/bigeye_level10_q_anchor.cpp"

if [[ ! -f "$DRIVER" ]]; then
  echo "ERROR: missing $DRIVER. Run from repo root after Level 10 patch."
  exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
cp "$DRIVER" "${DRIVER}.before_remove_stale_geometry.${STAMP}"

python3 - <<'PY'
from pathlib import Path
import re

p = Path("examples/NMFS/pifsc_bigeye_tuna/level10_q_anchor/quadra/bigeye_level10_q_anchor.cpp")
s = p.read_text()

# Remove any stale Level 9/10 fixed-effect geometry include.
s = re.sub(r'#include "\.\./diagnostics/bigeye_level[0-9]+_fixed_effect_geometry\.hpp"\n', '', s)

# Remove stale geometry report calls, including multi-line string literal args.
s = re.sub(
    r'\n\s*pifsc_bigeye_tuna::write_level[0-9]+_fixed_effect_geometry_report\(\s*\n'
    r'(?:.|\n)*?'
    r'\n\s*objective,\s*fit\);\s*\n',
    '\n',
    s,
    count=1,
)

s = re.sub(
    r'\n\s*pifsc_bigeye_tuna::write_level[0-9]+_fixed_effect_geometry_report\(\s*\n'
    r'(?:.|\n)*?'
    r'\n\s*objective,\s*params,\s*fit,\s*opts\);\s*\n',
    '\n',
    s,
    count=1,
)

# Remove stale cout lines for geometry outputs.
s = re.sub(r'\n\s*std::cout << "wrote:\s+.*fixed_effect_geometry_report\.(txt|csv)\\n";', '', s)

p.write_text(s)
PY

cat > inspect_bigeye_level10_remove_stale_geometry.sh <<'EOF_INSPECT'
#!/usr/bin/env bash
set -euo pipefail

DRIVER="examples/NMFS/pifsc_bigeye_tuna/level10_q_anchor/quadra/bigeye_level10_q_anchor.cpp"

echo "== Remaining stale geometry references? =="
grep -n "write_level9_fixed_effect_geometry_report\\|write_level10_fixed_effect_geometry_report\\|fixed_effect_geometry" "$DRIVER" || true

echo
echo "== Nearby report calls =="
grep -n "write_bigeye_report_suite\\|write_recruitment_diagnostics\\|objective_consistency" "$DRIVER" || true
EOF_INSPECT

chmod +x inspect_bigeye_level10_remove_stale_geometry.sh

cat > run_bigeye_level10_q_anchor_after_geometry_cleanup.sh <<'EOF_RUN'
#!/usr/bin/env bash
set -euo pipefail

./inspect_bigeye_level10_remove_stale_geometry.sh

echo
echo "== O3 build Bigeye Level 10 q anchor after geometry cleanup =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level10_q_anchor/quadra/bigeye_level10_q_anchor.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level10_q_anchor/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level10_q_anchor_check

echo
echo "== Run Bigeye Level 10 q anchor =="
./build/examples/pifsc_bigeye_level10_q_anchor_check

echo
echo "== Level 10 fit summary =="
cat examples/NMFS/pifsc_bigeye_tuna/level10_q_anchor/outputs/bigeye_level10_fit_summary.csv

echo
echo "== Level 10 components =="
cat examples/NMFS/pifsc_bigeye_tuna/level10_q_anchor/outputs/bigeye_level10_objective_components.csv

echo
echo "== Level 10 recruitment diagnostics =="
sed -n '1,150p' examples/NMFS/pifsc_bigeye_tuna/level10_q_anchor/outputs/bigeye_level10_recruitment_diagnostics.txt
EOF_RUN

chmod +x run_bigeye_level10_q_anchor_after_geometry_cleanup.sh

echo "Removed stale Level 9 geometry call from Level 10 driver."
echo
echo "Run:"
echo "  ./inspect_bigeye_level10_remove_stale_geometry.sh"
echo "  ./run_bigeye_level10_q_anchor_after_geometry_cleanup.sh"
