#!/usr/bin/env bash
set -euo pipefail

./inspect_bigeye_level2_fleet_specific_composition.sh

echo
echo "== O3 build Bigeye Level 2 fleet-specific composition =="
mkdir -p build/examples
c++ -std=c++17 -O3   -I.   -Iexternal/eigen   examples/NMFS/pifsc_bigeye_tuna/level2_fleet_specific_composition/quadra/bigeye_level2_fleet_specific_composition.cpp   examples/NMFS/pifsc_bigeye_tuna/level2_fleet_specific_composition/quadra/bigeye_adgraph_global.cpp   -o build/examples/pifsc_bigeye_level2_fleet_specific_composition_check

echo
echo "== Run Bigeye Level 2 fleet-specific composition =="
./build/examples/pifsc_bigeye_level2_fleet_specific_composition_check

echo
echo "== Level 2 fit summary =="
cat examples/NMFS/pifsc_bigeye_tuna/level2_fleet_specific_composition/outputs/bigeye_level2_fit_summary.csv

echo
echo "== Level 2 fixed-effect geometry preview =="
sed -n '1,180p'   examples/NMFS/pifsc_bigeye_tuna/level2_fleet_specific_composition/outputs/bigeye_level2_fixed_effect_geometry_report.txt
