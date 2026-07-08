#!/usr/bin/env bash
set -euo pipefail

./inspect_bigeye_level1b_fleet_specific_q.sh

echo
echo "== O3 build Bigeye Level 1B fleet-specific q =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level1b_fleet_specific_q/quadra/bigeye_level1b_fleet_specific_q.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level1b_fleet_specific_q/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level1b_fleet_specific_q_check

echo
echo "== Run Bigeye Level 1B fleet-specific q =="
./build/examples/pifsc_bigeye_level1b_fleet_specific_q_check

echo
echo "== Level 1B fit summary =="
cat examples/NMFS/pifsc_bigeye_tuna/level1b_fleet_specific_q/outputs/bigeye_level1b_fit_summary.csv

echo
echo "== Level 1B fixed-effect geometry preview =="
sed -n '1,180p' \
  examples/NMFS/pifsc_bigeye_tuna/level1b_fleet_specific_q/outputs/bigeye_level1b_fixed_effect_geometry_report.txt
