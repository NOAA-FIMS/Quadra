#!/usr/bin/env bash
set -euo pipefail

./inspect_bigeye_level6_purse_seine_age_selectivity.sh

echo
echo "== O3 build Bigeye Level 6 purse-seine age selectivity =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity/quadra/bigeye_level6_purse_seine_age_selectivity.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level6_purse_seine_age_selectivity_check

echo
echo "== Run Bigeye Level 6 purse-seine age selectivity =="
./build/examples/pifsc_bigeye_level6_purse_seine_age_selectivity_check

echo
echo "== Level 6 fit summary =="
cat examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity/outputs/bigeye_level6_fit_summary.csv

echo
echo "== Level 6 fixed-effect geometry preview =="
sed -n '1,240p' \
  examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity/outputs/bigeye_level6_fixed_effect_geometry_report.txt

echo
echo "== Level 6 wiggle diagnostics preview =="
sed -n '1,180p' \
  examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity/outputs/bigeye_level6_fixed_effect_wiggle_diagnostics.txt
