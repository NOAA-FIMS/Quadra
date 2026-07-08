#!/usr/bin/env bash
set -euo pipefail

./inspect_bigeye_level13_initial_numbers.sh

echo
echo "== O3 build Bigeye Level 13 estimated initial numbers =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level13_estimated_initial_numbers/quadra/bigeye_level13_estimated_initial_numbers.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level13_estimated_initial_numbers/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level13_initial_numbers_check

echo
echo "== Run Bigeye Level 13 estimated initial numbers =="
./build/examples/pifsc_bigeye_level13_initial_numbers_check

echo
echo "== Level 13 fit summary =="
cat examples/NMFS/pifsc_bigeye_tuna/level13_estimated_initial_numbers/outputs/bigeye_level13_fit_summary.csv

echo
echo "== Level 13 recruitment diagnostics =="
sed -n '1,150p' examples/NMFS/pifsc_bigeye_tuna/level13_estimated_initial_numbers/outputs/bigeye_level13_recruitment_diagnostics.txt

echo
echo "== Level 12 vs Level 13 recruitment comparison =="
python3 - <<'PY'
import csv
from pathlib import Path

base = Path("examples/NMFS/pifsc_bigeye_tuna")
paths = {
    "level12": base / "level12_tuna_life_history/outputs/bigeye_level12_recruitment_diagnostics.csv",
    "level13": base / "level13_estimated_initial_numbers/outputs/bigeye_level13_recruitment_diagnostics.csv",
}
metrics = ["sd", "lag1_correlation", "roughness", "total_prior_nll", "max_abs"]
print("level," + ",".join(metrics))
for level, p in paths.items():
    vals = {}
    if p.exists():
        with p.open() as f:
            for r in csv.DictReader(f):
                if r.get("section") == "summary" and r.get("metric") in metrics:
                    vals[r["metric"]] = r["value"]
    print(level + "," + ",".join(vals.get(m, "") for m in metrics))
PY
