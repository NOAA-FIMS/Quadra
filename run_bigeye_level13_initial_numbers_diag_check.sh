#!/usr/bin/env bash
set -euo pipefail

./inspect_bigeye_level13_initial_numbers_diag.sh

echo
echo "== O3 build Bigeye Level 13 initial numbers diagnostics =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level13_estimated_initial_numbers/quadra/bigeye_level13_estimated_initial_numbers.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level13_estimated_initial_numbers/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level13_initial_numbers_diag_check

echo
echo "== Run Bigeye Level 13 initial numbers diagnostics =="
./build/examples/pifsc_bigeye_level13_initial_numbers_diag_check

echo
echo "== Initial numbers diagnostics preview =="
sed -n '1,120p' \
  examples/NMFS/pifsc_bigeye_tuna/level13_estimated_initial_numbers/outputs/bigeye_level13_initial_numbers_diagnostics.txt

echo
echo "== Compact old-fish summary =="
python3 - <<'PY'
import csv
from pathlib import Path

p = Path("examples/NMFS/pifsc_bigeye_tuna/level13_estimated_initial_numbers/outputs/bigeye_level13_initial_numbers_diagnostics.csv")
vals = {}
rows = {}
with p.open() as f:
    for r in csv.DictReader(f):
        if r["section"] == "summary":
            vals[r["metric"]] = float(r["value"])
        elif r["section"] == "age":
            rows.setdefault(r["target"], {})[r["metric"]] = float(r["value"])

print("summary_metric,value")
for k in [
    "equilibrium_plus_n_share",
    "fitted_plus_n_share",
    "equilibrium_plus_bio_share",
    "fitted_plus_bio_share",
    "fitted_over_equilibrium_n",
    "fitted_over_equilibrium_bio",
    "fitted_over_equilibrium_ssb",
]:
    print(f"{k},{vals.get(k, float('nan'))}")

print()
print("age,multiplier,equilibrium_biomass_share,fitted_biomass_share")
for age in ["age_8", "age_9", "age_10"]:
    r = rows.get(age, {})
    print(f"{age},{r.get('multiplier', float('nan'))},{r.get('equilibrium_biomass_share', float('nan'))},{r.get('fitted_biomass_share', float('nan'))}")
PY
