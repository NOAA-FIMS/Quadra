#!/usr/bin/env bash
set -euo pipefail

./inspect_bigeye_level7_dual_age_selectivity.sh

echo
echo "== O3 build Bigeye Level 7 dual age selectivity =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level7_dual_age_selectivity/quadra/bigeye_level7_dual_age_selectivity.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level7_dual_age_selectivity/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level7_dual_age_selectivity_check

echo
echo "== Run Bigeye Level 7 dual age selectivity =="
./build/examples/pifsc_bigeye_level7_dual_age_selectivity_check

echo
echo "== Level 7 fit summary =="
cat examples/NMFS/pifsc_bigeye_tuna/level7_dual_age_selectivity/outputs/bigeye_level7_fit_summary.csv

echo
echo "== Level 7 objective components =="
cat examples/NMFS/pifsc_bigeye_tuna/level7_dual_age_selectivity/outputs/bigeye_level7_objective_components.csv

echo
echo "== Level 7 recruitment diagnostics preview =="
sed -n '1,180p' \
  examples/NMFS/pifsc_bigeye_tuna/level7_dual_age_selectivity/outputs/bigeye_level7_recruitment_diagnostics.txt

echo
echo "== Level 6 vs Level 7 recruitment comparison =="
python3 - <<'PY'
import csv
from pathlib import Path

base = Path("examples/NMFS/pifsc_bigeye_tuna")

paths = {
    "level6": base / "level6_purse_seine_age_selectivity/outputs/bigeye_level6_recruitment_diagnostics.csv",
    "level7": base / "level7_dual_age_selectivity/outputs/bigeye_level7_recruitment_diagnostics.csv",
}

metrics = ["sd", "lag1_correlation", "roughness", "total_prior_nll", "max_abs"]

print("level," + ",".join(metrics))
for level, p in paths.items():
    vals = {}
    if not p.exists():
        continue
    with p.open() as f:
        for r in csv.DictReader(f):
            if r["section"] == "summary" and r["metric"] in metrics:
                vals[r["metric"]] = r["value"]
    print(level + "," + ",".join(vals.get(m, "") for m in metrics))
PY

echo
echo "== Recommended commit status =="
git status --short examples/NMFS/pifsc_bigeye_tuna/level7_dual_age_selectivity examples/NMFS/pifsc_bigeye_tuna/workflow/scientific_reasoning_log.md
