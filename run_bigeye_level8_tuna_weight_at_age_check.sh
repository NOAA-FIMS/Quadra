#!/usr/bin/env bash
set -euo pipefail

./inspect_bigeye_level8_tuna_weight_at_age.sh

echo
echo "== O3 build Bigeye Level 8 tuna weight-at-age =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level8_tuna_weight_at_age/quadra/bigeye_level8_tuna_weight_at_age.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level8_tuna_weight_at_age/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level8_tuna_weight_at_age_check

echo
echo "== Run Bigeye Level 8 tuna weight-at-age =="
./build/examples/pifsc_bigeye_level8_tuna_weight_at_age_check

echo
echo "== Level 8 fit summary =="
cat examples/NMFS/pifsc_bigeye_tuna/level8_tuna_weight_at_age/outputs/bigeye_level8_fit_summary.csv

echo
echo "== Level 8 objective components =="
cat examples/NMFS/pifsc_bigeye_tuna/level8_tuna_weight_at_age/outputs/bigeye_level8_objective_components.csv || true

echo
echo "== Level 8 recruitment diagnostics =="
sed -n '1,170p' \
  examples/NMFS/pifsc_bigeye_tuna/level8_tuna_weight_at_age/outputs/bigeye_level8_recruitment_diagnostics.txt

echo
echo "== Level 6 vs Level 8 recruitment comparison =="
python3 - <<'PY'
import csv
from pathlib import Path

base = Path("examples/NMFS/pifsc_bigeye_tuna")
paths = {
    "level6": base / "level6_purse_seine_age_selectivity/outputs/bigeye_level6_recruitment_diagnostics.csv",
    "level8": base / "level8_tuna_weight_at_age/outputs/bigeye_level8_recruitment_diagnostics.csv",
}
metrics = ["sd", "lag1_correlation", "roughness", "total_prior_nll", "max_abs"]
print("level," + ",".join(metrics))
for level, p in paths.items():
    vals = {}
    if p.exists():
        with p.open() as f:
            for r in csv.DictReader(f):
                if r["section"] == "summary" and r["metric"] in metrics:
                    vals[r["metric"]] = r["value"]
    print(level + "," + ",".join(vals.get(m, "") for m in metrics))
PY
