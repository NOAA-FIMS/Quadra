#!/usr/bin/env bash
set -euo pipefail

./inspect_bigeye_level6_safe_diagnostics.sh

echo
echo "== O3 build Bigeye Level 6 safe diagnostics =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity/quadra/bigeye_level6_purse_seine_age_selectivity.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level6_safe_diagnostics_check

echo
echo "== Run Bigeye Level 6 safe diagnostics =="
./build/examples/pifsc_bigeye_level6_safe_diagnostics_check

echo
echo "== Safe wiggle diagnostics preview =="
sed -n '1,180p' \
  examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity/outputs/bigeye_level6_safe_fixed_effect_wiggle_diagnostics.txt
