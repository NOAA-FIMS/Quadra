#!/usr/bin/env bash
set -euo pipefail

./inspect_bigeye_level12_tuna_life_history.sh

echo
echo "== O3 build Bigeye Level 12 tuna life history =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level12_tuna_life_history/quadra/bigeye_level12_tuna_life_history.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level12_tuna_life_history/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level12_tuna_life_history_check

echo
echo "== Run Bigeye Level 12 tuna life history =="
./build/examples/pifsc_bigeye_level12_tuna_life_history_check

echo
echo "== Level 12 fit summary =="
cat examples/NMFS/pifsc_bigeye_tuna/level12_tuna_life_history/outputs/bigeye_level12_fit_summary.csv

echo
echo "== Level 12 components =="
cat examples/NMFS/pifsc_bigeye_tuna/level12_tuna_life_history/outputs/bigeye_level12_objective_components.csv

echo
echo "== Level 12 recruitment diagnostics =="
sed -n '1,150p' examples/NMFS/pifsc_bigeye_tuna/level12_tuna_life_history/outputs/bigeye_level12_recruitment_diagnostics.txt

echo
echo "== Level 11 vs Level 12 recruitment comparison =="
python3 - <<'PY'
import csv
from pathlib import Path

base = Path("examples/NMFS/pifsc_bigeye_tuna")
paths = {
    "level11": base / "level11_fixed_m_q_anchor/outputs/bigeye_level11_recruitment_diagnostics.csv",
    "level12": base / "level12_tuna_life_history/outputs/bigeye_level12_recruitment_diagnostics.csv",
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
