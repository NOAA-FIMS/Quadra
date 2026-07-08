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
