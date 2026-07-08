#!/usr/bin/env bash
set -euo pipefail
./inspect_bigeye_level16_purse_seine_age_selectivity.sh
echo
echo "== O3 build Bigeye Level 16 purse-seine age selectivity =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level16_purse_seine_age_selectivity/quadra/bigeye_level16_purse_seine_age_selectivity.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level16_purse_seine_age_selectivity/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level16_purse_seine_age_selectivity_check
echo
echo "== Run Bigeye Level 16 purse-seine age selectivity =="
./build/examples/pifsc_bigeye_level16_purse_seine_age_selectivity_check
echo
echo "== Level 16 fit summary =="
cat examples/NMFS/pifsc_bigeye_tuna/level16_purse_seine_age_selectivity/outputs/bigeye_level16_fit_summary.csv
echo
echo "== Purse-seine selectivity estimates =="
grep 'sel_purse_seine_age_' examples/NMFS/pifsc_bigeye_tuna/level16_purse_seine_age_selectivity/outputs/bigeye_level16_fit_summary.csv || true
echo
echo "== Level 16 age-comp residual summary =="
grep -n -A20 'Fleet summary' examples/NMFS/pifsc_bigeye_tuna/level16_purse_seine_age_selectivity/outputs/bigeye_level16_age_comp_residual_diagnostics.txt || true
