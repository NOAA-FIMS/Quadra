#!/usr/bin/env bash
set -euo pipefail

echo "== O3 build Bigeye Level 19 parameter sanity diagnostics =="
mkdir -p build/examples
c++ -std=c++17 -O3 -flto \
  -I. -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level19_flexible_longline_selectivity/quadra/bigeye_level19_flexible_longline_selectivity.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level19_flexible_longline_selectivity/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level19_parameter_sanity_check

echo
echo "== Run Bigeye Level 19 parameter sanity diagnostics =="
./build/examples/pifsc_bigeye_level19_parameter_sanity_check

echo
echo "== Level 19 fit summary key rows =="
grep -E "objective|grad_norm|converged|sel_longline_age_|init_number_multiplier_age_" \
  examples/NMFS/pifsc_bigeye_tuna/level19_flexible_longline_selectivity/outputs/bigeye_level19_fit_summary.csv \
  | head -90

echo
echo "== Level 19 prior penalty blocks =="
grep -n "Prior penalty by block" -A10 \
  examples/NMFS/pifsc_bigeye_tuna/level19_flexible_longline_selectivity/outputs/bigeye_level19_parameter_sanity_diagnostics.txt

echo
echo "== Level 19 initial-number stress =="
grep -n "Initial-number stress summary" -A10 \
  examples/NMFS/pifsc_bigeye_tuna/level19_flexible_longline_selectivity/outputs/bigeye_level19_parameter_sanity_diagnostics.txt

echo
echo "== Level 19 age-specific parameter sanity =="
grep -n "Age-specific parameters" -A16 \
  examples/NMFS/pifsc_bigeye_tuna/level19_flexible_longline_selectivity/outputs/bigeye_level19_parameter_sanity_diagnostics.txt

echo
echo "== Level 18 vs Level 19 objective/residual comparison =="
python3 - <<'PY'
from pathlib import Path
import csv

def read_summary(path):
    out = {}
    if not Path(path).exists():
        return out
    with open(path) as f:
        for row in csv.reader(f):
            if len(row) >= 2:
                out[row[0]] = row[1]
    return out

def fleet_summary(path):
    out = {}
    if not Path(path).exists():
        return out
    lines = Path(path).read_text().splitlines()
    in_block = False
    for line in lines:
        if line.startswith("fleet,mean_abs_residual"):
            in_block = True
            continue
        if in_block:
            if not line.strip():
                break
            parts = line.split(",")
            if len(parts) >= 4:
                out[parts[0]] = {"mean_abs": parts[1], "max_abs": parts[2], "n": parts[3]}
    return out

s18 = read_summary("examples/NMFS/pifsc_bigeye_tuna/level18_longline_selectivity_diagnostic/outputs/bigeye_level18_fit_summary.csv")
s19 = read_summary("examples/NMFS/pifsc_bigeye_tuna/level19_flexible_longline_selectivity/outputs/bigeye_level19_fit_summary.csv")
f18 = fleet_summary("examples/NMFS/pifsc_bigeye_tuna/level18_longline_selectivity_diagnostic/outputs/bigeye_level18_age_comp_residual_diagnostics.txt")
f19 = fleet_summary("examples/NMFS/pifsc_bigeye_tuna/level19_flexible_longline_selectivity/outputs/bigeye_level19_age_comp_residual_diagnostics.txt")

print("level,objective,grad_norm,longline_mean_abs,longline_max_abs,purse_seine_mean_abs,purse_seine_max_abs")
print(f"level18,{s18.get('objective','')},{s18.get('grad_norm','')},{f18.get('longline',{}).get('mean_abs','')},{f18.get('longline',{}).get('max_abs','')},{f18.get('purse_seine',{}).get('mean_abs','')},{f18.get('purse_seine',{}).get('max_abs','')}")
print(f"level19,{s19.get('objective','')},{s19.get('grad_norm','')},{f19.get('longline',{}).get('mean_abs','')},{f19.get('longline',{}).get('max_abs','')},{f19.get('purse_seine',{}).get('mean_abs','')},{f19.get('purse_seine',{}).get('max_abs','')}")
PY
