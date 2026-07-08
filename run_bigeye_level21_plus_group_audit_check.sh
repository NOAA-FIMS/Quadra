#!/usr/bin/env bash
set -euo pipefail

SRC="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/quadra/bigeye_level21_age_based_natural_mortality_diagnostic.cpp"
ADG="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/quadra/bigeye_adgraph_global.cpp"
BIN="/tmp/bigeye_level21_plus_group_audit"

echo "== O3 build Bigeye Level 21 plus-group audit =="
clang++ -std=c++17 -O3 -Iexternal/eigen -I. "$SRC" "$ADG" -o "$BIN"

echo
echo "== Run Bigeye Level 21 plus-group audit =="
"$BIN"

OUT="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_plus_group_audit.txt"

echo
echo "== Plus-group initialization =="
grep -A12 "Initialization" "$OUT"

echo
echo "== Initial state tail rows =="
grep -A14 "Initial state by age" "$OUT" | tail -6

echo
echo "== Annual plus-group dynamics preview =="
grep -A25 "Annual plus-group dynamics" "$OUT" | head -28

echo
echo "== Longline tail comp top residuals from CSV =="
python3 - <<'PY'
from pathlib import Path
import csv
p=Path("examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_plus_group_audit.csv")
rows=[]
with p.open() as f:
    r=csv.DictReader(f)
    for row in r:
        if row["section"]=="tail_comp":
            try:
                row["abs_resid"]=abs(float(row["residual"]))
                rows.append(row)
            except Exception:
                pass
rows.sort(key=lambda x:x["abs_resid"], reverse=True)
print("year,age,n,ll_sel,pred_ll_comp,obs_ll_comp,residual,abs_residual")
for row in rows[:20]:
    print(",".join([row["year"], row["age"], row["n"], row["ll_sel"], row["pred_ll_comp"], row["obs_ll_comp"], row["residual"], str(row["abs_resid"])]))
PY
