#!/usr/bin/env bash
set -euo pipefail

./inspect_bigeye_level0_clone.sh

echo
echo "== O3 build Bigeye Level 0 clone =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level0_single_region/quadra/bigeye_level0.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level0_single_region/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level0_check

echo
echo "== Run Bigeye Level 0 clone =="
./build/examples/pifsc_bigeye_level0_check

echo
echo "== Bigeye Level 0 outputs =="
find examples/NMFS/pifsc_bigeye_tuna/level0_single_region/outputs \
  -maxdepth 1 -type f | sort
