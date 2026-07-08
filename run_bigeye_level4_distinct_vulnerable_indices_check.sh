#!/usr/bin/env bash
set -euo pipefail

./inspect_bigeye_level4_distinct_vulnerable_indices.sh

echo
echo "== O3 build Bigeye Level 4 distinct vulnerable indices =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level4_distinct_vulnerable_indices/quadra/bigeye_level4_distinct_vulnerable_indices.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level4_distinct_vulnerable_indices/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level4_distinct_vulnerable_indices_check

echo
echo "== Run Bigeye Level 4 distinct vulnerable indices =="
./build/examples/pifsc_bigeye_level4_distinct_vulnerable_indices_check

echo
echo "== Level 4 fit summary =="
cat examples/NMFS/pifsc_bigeye_tuna/level4_distinct_vulnerable_indices/outputs/bigeye_level4_fit_summary.csv

echo
echo "== Level 4 fixed-effect geometry preview =="
sed -n '1,260p' \
  examples/NMFS/pifsc_bigeye_tuna/level4_distinct_vulnerable_indices/outputs/bigeye_level4_fixed_effect_geometry_report.txt
