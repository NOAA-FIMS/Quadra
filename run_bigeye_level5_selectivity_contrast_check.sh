#!/usr/bin/env bash
set -euo pipefail

./inspect_bigeye_level5_selectivity_contrast.sh

echo
echo "== O3 build Bigeye Level 5 selectivity contrast =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level5_selectivity_contrast/quadra/bigeye_level5_selectivity_contrast.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level5_selectivity_contrast/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level5_selectivity_contrast_check

echo
echo "== Run Bigeye Level 5 selectivity contrast =="
./build/examples/pifsc_bigeye_level5_selectivity_contrast_check

echo
echo "== Level 5 fit summary =="
cat examples/NMFS/pifsc_bigeye_tuna/level5_selectivity_contrast/outputs/bigeye_level5_fit_summary.csv

echo
echo "== Level 5 fixed-effect geometry preview =="
sed -n '1,260p' \
  examples/NMFS/pifsc_bigeye_tuna/level5_selectivity_contrast/outputs/bigeye_level5_fixed_effect_geometry_report.txt

echo
echo "== Level 5 wiggle diagnostics preview =="
sed -n '1,180p' \
  examples/NMFS/pifsc_bigeye_tuna/level5_selectivity_contrast/outputs/bigeye_level5_fixed_effect_wiggle_diagnostics.txt
