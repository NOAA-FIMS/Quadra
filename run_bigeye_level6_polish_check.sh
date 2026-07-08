#!/usr/bin/env bash
set -euo pipefail

./inspect_bigeye_level6_polish.sh

echo
echo "== O3 build Bigeye Level 6 polished diagnostics =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity/quadra/bigeye_level6_purse_seine_age_selectivity.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level6_polish_check

echo
echo "== Run Bigeye Level 6 polished diagnostics =="
./build/examples/pifsc_bigeye_level6_polish_check

echo
echo "== Safe geometry summary preview =="
sed -n '1,120p' \
  examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity/outputs/bigeye_level6_safe_fixed_effect_geometry_summary.txt

echo
echo "== Safe wiggle diagnostics preview =="
sed -n '1,120p' \
  examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity/outputs/bigeye_level6_safe_fixed_effect_wiggle_diagnostics.txt

echo
echo "== Recommended commit status =="
git status --short examples/NMFS/pifsc_bigeye_tuna/level6_purse_seine_age_selectivity examples/NMFS/pifsc_bigeye_tuna/workflow/scientific_reasoning_log.md
