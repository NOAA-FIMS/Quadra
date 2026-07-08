#!/usr/bin/env bash
set -euo pipefail

./inspect_bigeye_level6_recruitment_diagnostics.sh

echo
echo "== O3 build Bigeye Level 6 recruitment diagnostics =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity/quadra/bigeye_level6_purse_seine_age_selectivity.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level6_recruitment_diagnostics_check

echo
echo "== Run Bigeye Level 6 recruitment diagnostics =="
./build/examples/pifsc_bigeye_level6_recruitment_diagnostics_check

echo
echo "== Recruitment diagnostics preview =="
sed -n '1,180p' \
  examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity/outputs/bigeye_level6_recruitment_diagnostics.txt

echo
echo "== Objective components =="
cat examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity/outputs/bigeye_level6_objective_components.csv

echo
echo "== Top recruitment deviations by abs value =="
python3 - <<'PY'
import csv

p = "examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity/outputs/bigeye_level6_recruitment_diagnostics.csv"

rows = []
with open(p) as f:
    for r in csv.DictReader(f):
        if r["section"] == "recruitment" and r["metric"] == "rec_dev":
            rows.append((abs(float(r["value"])), r["target"], float(r["value"])))

rows.sort(reverse=True)
print("year,rec_dev,abs_rec_dev")
for abs_v, year, v in rows[:10]:
    print(f"{year},{v:.6g},{abs_v:.6g}")
PY
