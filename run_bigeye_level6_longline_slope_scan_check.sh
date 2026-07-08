#!/usr/bin/env bash
set -euo pipefail

./inspect_bigeye_level6_longline_slope_scan.sh

echo
echo "== O3 build Bigeye Level 6 longline slope scan =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity/quadra/bigeye_level6_purse_seine_age_selectivity.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level6_longline_slope_scan_check

echo
echo "== Run Bigeye Level 6 longline slope scan =="
./build/examples/pifsc_bigeye_level6_longline_slope_scan_check

echo
echo "== Longline slope geometry scan preview =="
sed -n '1,180p' \
  examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity/outputs/bigeye_level6_longline_slope_geometry_scan.txt

echo
echo "== Scan compact table =="
python3 - <<'PY'
import csv

p = "examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity/outputs/bigeye_level6_longline_slope_geometry_scan.csv"

with open(p) as f:
    rows = list(csv.DictReader(f))

print("multiplier,sel_slope,delta_objective,status,message")
for r in rows:
    print(",".join([
        f'{float(r["multiplier"]):.3f}',
        f'{float(r["sel_slope_longline"]):.6g}',
        r["delta_objective"],
        r["status"],
        r["message"][:80],
    ]))
PY
