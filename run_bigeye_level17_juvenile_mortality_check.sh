#!/usr/bin/env bash
set -euo pipefail

./inspect_bigeye_level17_juvenile_mortality.sh

echo
echo "== O3 build Bigeye Level 17 juvenile mortality diagnostic =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level17_juvenile_mortality_diagnostic/quadra/bigeye_level17_juvenile_mortality_diagnostic.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level17_juvenile_mortality_diagnostic/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level17_juvenile_mortality_check

echo
echo "== Run Bigeye Level 17 juvenile mortality diagnostic =="
./build/examples/pifsc_bigeye_level17_juvenile_mortality_check

echo
echo "== Level 17 fit summary =="
cat examples/NMFS/pifsc_bigeye_tuna/level17_juvenile_mortality_diagnostic/outputs/bigeye_level17_fit_summary.csv

echo
echo "== Level 17 age-comp residual summary =="
grep -n -A20 "Fleet summary" examples/NMFS/pifsc_bigeye_tuna/level17_juvenile_mortality_diagnostic/outputs/bigeye_level17_age_comp_residual_diagnostics.txt || true

echo
echo "== Level 16 vs Level 17 compact comparison =="
python3 - <<'PY'
import csv
from pathlib import Path
base=Path("examples/NMFS/pifsc_bigeye_tuna")
def fit(level,prefix):
    d={}
    p=base/level/"outputs"/f"{prefix}_fit_summary.csv"
    if p.exists():
        for r in csv.DictReader(p.open()): d[r["field"]]=r["value"]
    return d
def rec(level,prefix):
    d={}
    p=base/level/"outputs"/f"{prefix}_recruitment_diagnostics.csv"
    if p.exists():
        for r in csv.DictReader(p.open()):
            if r.get("section")=="summary": d[r["metric"]]=r["value"]
    return d
def fleet(level,prefix):
    d={}
    p=base/level/"outputs"/f"{prefix}_age_comp_residual_diagnostics.csv"
    if p.exists():
        for r in csv.DictReader(p.open()):
            if r["section"]=="fleet_summary": d[r["fleet"]]=r["residual"]
    return d
items=[("level16","level16_purse_seine_age_selectivity","bigeye_level16"),
       ("level17","level17_juvenile_mortality_diagnostic","bigeye_level17")]
print("level,objective,grad_norm,juvenile_m_multiplier,juvenile_m,rec_sd,rec_lag1,rec_prior_nll,purse_seine_mean_abs,longline_mean_abs")
for label,level,prefix in items:
    fi=fit(level,prefix); re=rec(level,prefix); fl=fleet(level,prefix)
    print(",".join([label,fi.get("objective",""),fi.get("grad_norm",""),
                    fi.get("juvenile_m_multiplier",""),fi.get("juvenile_m_age1_2",""),
                    re.get("sd",""),re.get("lag1_correlation",""),re.get("total_prior_nll",""),
                    fl.get("purse_seine",""),fl.get("longline","")]))
PY
