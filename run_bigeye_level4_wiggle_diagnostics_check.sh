#!/usr/bin/env bash
set -euo pipefail

./inspect_bigeye_level4_wiggle_diagnostics.sh

echo
echo "== O3 build Bigeye Level 4 wiggle diagnostics =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level4_distinct_vulnerable_indices/quadra/bigeye_level4_distinct_vulnerable_indices.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level4_distinct_vulnerable_indices/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level4_wiggle_diagnostics_check

echo
echo "== Run Bigeye Level 4 wiggle diagnostics =="
./build/examples/pifsc_bigeye_level4_wiggle_diagnostics_check

echo
echo "== Wiggle diagnostics text preview =="
sed -n '1,180p' \
  examples/NMFS/pifsc_bigeye_tuna/level4_distinct_vulnerable_indices/outputs/bigeye_level4_fixed_effect_wiggle_diagnostics.txt
