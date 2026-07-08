#!/usr/bin/env bash
set -euo pipefail

L9="examples/NMFS/pifsc_bigeye_tuna/level9_estimated_natural_mortality"
REPORT="$L9/reports/bigeye_fit_reports.hpp"
SUITE="$L9/reports/bigeye_report_suite.hpp"

if [[ ! -f "$REPORT" || ! -f "$SUITE" ]]; then
  echo "ERROR: missing Level 9 report files. Run from repo root."
  exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
cp "$REPORT" "${REPORT}.before_report_compile_fix2.${STAMP}"
cp "$SUITE" "${SUITE}.before_report_compile_fix2.${STAMP}"

python3 - <<'PY'
from pathlib import Path

p = Path("examples/NMFS/pifsc_bigeye_tuna/level9_estimated_natural_mortality/reports/bigeye_fit_reports.hpp")
s = p.read_text()

# OptResult does not expose a raw joint member in this branch. The direct
# objective consistency check already writes the true joint objective at
# fit.par + fit.u_hat, so the report summary should not depend on a missing
# OptResult field.
s = s.replace('  out << "joint_objective," << fit.joint << "\\n";\n', '')
s = s.replace('  out << "joint_objective," << fit.joint << ",from optimizer\\n";\n', '')

# Defensive cleanup in case a previous patch used a different spelling.
s = s.replace('  out << "joint_objective," << fit.joint_value << "\\n";\n', '')
s = s.replace('  out << "joint_objective," << fit.joint_value << ",from optimizer\\n";\n', '')

p.write_text(s)
PY

# The Level 9 report implementation now lives entirely in bigeye_fit_reports.hpp.
# Keep this wrapper as a non-defining compatibility include to avoid duplicate
# default_bigeye_report_paths/write_bigeye_report_suite definitions.
cat > "$SUITE" <<'EOF_SUITE'
#pragma once

#include "bigeye_fit_reports.hpp"
EOF_SUITE

cat > inspect_bigeye_level9_report_compile_fix2.sh <<'EOF_INSPECT'
#!/usr/bin/env bash
set -euo pipefail

L9="examples/NMFS/pifsc_bigeye_tuna/level9_estimated_natural_mortality"
REPORT="$L9/reports/bigeye_fit_reports.hpp"
SUITE="$L9/reports/bigeye_report_suite.hpp"

echo "== Missing OptResult member references? =="
grep -n "fit\\.joint\\|fit\\.joint_value" "$REPORT" || true

echo
echo "== Suite wrapper =="
cat "$SUITE"

echo
echo "== Single definitions =="
grep -R "inline .*default_bigeye_report_paths\\|inline void write_bigeye_report_suite" -n "$L9/reports"
EOF_INSPECT

chmod +x inspect_bigeye_level9_report_compile_fix2.sh

cat > run_bigeye_level9_report_compile_fix2_check.sh <<'EOF_RUN'
#!/usr/bin/env bash
set -euo pipefail

./inspect_bigeye_level9_report_compile_fix2.sh

echo
echo "== O3 build Bigeye Level 9 report compile fix 2 =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level9_estimated_natural_mortality/quadra/bigeye_level9_estimated_natural_mortality.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level9_estimated_natural_mortality/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level9_report_compile_fix2_check

echo
echo "== Run Bigeye Level 9 report compile fix 2 =="
./build/examples/pifsc_bigeye_level9_report_compile_fix2_check

echo
echo "== Objective consistency check =="
cat examples/NMFS/pifsc_bigeye_tuna/level9_estimated_natural_mortality/outputs/bigeye_level9_objective_consistency_check.csv

echo
echo "== Patched objective components =="
cat examples/NMFS/pifsc_bigeye_tuna/level9_estimated_natural_mortality/outputs/bigeye_level9_objective_components.csv

echo
echo "== Consistency comparison =="
python3 - <<'PY'
import csv
from pathlib import Path

cons = Path("examples/NMFS/pifsc_bigeye_tuna/level9_estimated_natural_mortality/outputs/bigeye_level9_objective_consistency_check.csv")
comp = Path("examples/NMFS/pifsc_bigeye_tuna/level9_estimated_natural_mortality/outputs/bigeye_level9_objective_components.csv")

vals = {}
with cons.open() as f:
    for r in csv.DictReader(f):
        vals[r["metric"]] = float(r["value"])

components = {}
with comp.open() as f:
    for r in csv.DictReader(f):
        components[r["component"]] = float(r["value"])

print(f'direct_joint_objective={vals["direct_joint_objective"]:.12f}')
print(f'component_joint_total={components["joint_total"]:.12f}')
print(f'difference={components["joint_total"] - vals["direct_joint_objective"]:.12g}')
PY
EOF_RUN

chmod +x run_bigeye_level9_report_compile_fix2_check.sh

echo "Patched Level 9 report compile errors, pass 2."
echo
echo "Run:"
echo "  ./inspect_bigeye_level9_report_compile_fix2.sh"
echo "  ./run_bigeye_level9_report_compile_fix2_check.sh"
