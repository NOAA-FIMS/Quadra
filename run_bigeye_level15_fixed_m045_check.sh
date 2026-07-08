#!/usr/bin/env bash
set -euo pipefail

./inspect_bigeye_level15_fixed_m045.sh

echo
echo "== O3 build Bigeye Level 15 fixed M=0.45 best diagnostic =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level15_fixed_m045_best_diagnostic/quadra/bigeye_level15_fixed_m045_best_diagnostic.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level15_fixed_m045_best_diagnostic/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level15_fixed_m045_check

echo
echo "== Run Bigeye Level 15 fixed M=0.45 best diagnostic =="
./build/examples/pifsc_bigeye_level15_fixed_m045_check

echo
echo "== Level 15 fit summary =="
cat examples/NMFS/pifsc_bigeye_tuna/level15_fixed_m045_best_diagnostic/outputs/bigeye_level15_fit_summary.csv

echo
echo "== Level 15 recruitment diagnostics =="
sed -n '1,120p' examples/NMFS/pifsc_bigeye_tuna/level15_fixed_m045_best_diagnostic/outputs/bigeye_level15_recruitment_diagnostics.txt

echo
echo "== Level 15 initial numbers diagnostics =="
sed -n '1,90p' examples/NMFS/pifsc_bigeye_tuna/level15_fixed_m045_best_diagnostic/outputs/bigeye_level15_initial_numbers_diagnostics.txt

echo
echo "== Level 15 age-comp residual summary =="
grep -n -A20 "Fleet summary" examples/NMFS/pifsc_bigeye_tuna/level15_fixed_m045_best_diagnostic/outputs/bigeye_level15_age_comp_residual_diagnostics.txt || true
grep -n -A25 "Worst years" examples/NMFS/pifsc_bigeye_tuna/level15_fixed_m045_best_diagnostic/outputs/bigeye_level15_age_comp_residual_diagnostics.txt || true

echo
echo "== Compact worst age-comp residuals =="
python3 - <<'PY'
import csv
from pathlib import Path

p = Path("examples/NMFS/pifsc_bigeye_tuna/level15_fixed_m045_best_diagnostic/outputs/bigeye_level15_age_comp_residual_diagnostics.csv")
rows = []
with p.open() as f:
    for r in csv.DictReader(f):
        if r["section"] == "age_residual":
            try:
                r["_abs"] = float(r["abs_residual"])
                rows.append(r)
            except Exception:
                pass

print("year,fleet,age,observed,predicted,residual,abs_residual")
for r in sorted(rows, key=lambda x: x["_abs"], reverse=True)[:20]:
    print(",".join([r["year"], r["fleet"], r["age"], r["observed"], r["predicted"], r["residual"], r["abs_residual"]]))
PY
