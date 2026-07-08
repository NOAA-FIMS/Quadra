#!/usr/bin/env bash
set -euo pipefail

./inspect_bigeye_level3_fleet_specific_selectivity_index.sh

echo
echo "== O3 build Bigeye Level 3 fleet-specific selectivity/index =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level3_fleet_specific_selectivity_index/quadra/bigeye_level3_fleet_specific_selectivity_index.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level3_fleet_specific_selectivity_index/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level3_fleet_specific_selectivity_index_check

echo
echo "== Run Bigeye Level 3 fleet-specific selectivity/index =="
./build/examples/pifsc_bigeye_level3_fleet_specific_selectivity_index_check

echo
echo "== Level 3 fit summary =="
cat examples/NMFS/pifsc_bigeye_tuna/level3_fleet_specific_selectivity_index/outputs/bigeye_level3_fit_summary.csv

echo
echo "== Level 3 fixed-effect geometry preview =="
sed -n '1,240p' \
  examples/NMFS/pifsc_bigeye_tuna/level3_fleet_specific_selectivity_index/outputs/bigeye_level3_fixed_effect_geometry_report.txt
