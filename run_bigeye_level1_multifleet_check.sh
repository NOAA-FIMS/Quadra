#!/usr/bin/env bash
set -euo pipefail

./inspect_bigeye_level1_multifleet.sh

echo
echo "== O3 build Bigeye Level 1 multi-fleet scaffold =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level1_multifleet/quadra/bigeye_level1_multifleet.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level1_multifleet/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level1_multifleet_check

echo
echo "== Run Bigeye Level 1 multi-fleet scaffold =="
./build/examples/pifsc_bigeye_level1_multifleet_check

echo
echo "== Bigeye Level 1 outputs =="
find examples/NMFS/pifsc_bigeye_tuna/level1_multifleet/outputs \
  -maxdepth 1 -type f | sort

echo
echo "== Level 1 fit summary =="
cat examples/NMFS/pifsc_bigeye_tuna/level1_multifleet/outputs/bigeye_level1_fit_summary.csv

echo
echo "== Level 1 Laplace structure preview =="
sed -n '1,80p' \
  examples/NMFS/pifsc_bigeye_tuna/level1_multifleet/outputs/bigeye_level1_laplace_structure_report.txt
