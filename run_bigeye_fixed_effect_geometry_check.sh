#!/usr/bin/env bash
set -euo pipefail

./inspect_bigeye_fixed_effect_geometry.sh

echo
echo "== O3 build Bigeye fixed-effect geometry =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level1_multifleet/quadra/bigeye_level1_multifleet.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level1_multifleet/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_fixed_effect_geometry_check

echo
echo "== Run Bigeye fixed-effect geometry check =="
./build/examples/pifsc_bigeye_fixed_effect_geometry_check

echo
echo "== Fixed-effect geometry text preview =="
sed -n '1,180p' \
  examples/NMFS/pifsc_bigeye_tuna/level1_multifleet/outputs/bigeye_level1_fixed_effect_geometry_report.txt

echo
echo "== Fixed-effect geometry CSV preview =="
sed -n '1,80p' \
  examples/NMFS/pifsc_bigeye_tuna/level1_multifleet/outputs/bigeye_level1_fixed_effect_geometry_report.csv
