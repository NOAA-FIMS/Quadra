#!/usr/bin/env bash
set -euo pipefail

echo "== O3 build Bigeye Level 18 longline selectivity diagnostic =="
mkdir -p build/examples
c++ -std=c++17 -O3 -flto \
  -I. -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level18_longline_selectivity_diagnostic/quadra/bigeye_level18_longline_selectivity_diagnostic.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level18_longline_selectivity_diagnostic/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level18_longline_selectivity_diagnostic_check

echo
echo "== Run Bigeye Level 18 longline selectivity diagnostic =="
./build/examples/pifsc_bigeye_level18_longline_selectivity_diagnostic_check

echo
echo "== Level 18 fit summary =="
cat examples/NMFS/pifsc_bigeye_tuna/level18_longline_selectivity_diagnostic/outputs/bigeye_level18_fit_summary.csv

echo
echo "== Objective-consistent age-comp residual summary =="
grep -n "Fleet summary" -A8 \
  examples/NMFS/pifsc_bigeye_tuna/level18_longline_selectivity_diagnostic/outputs/bigeye_level18_age_comp_residual_diagnostics.txt

echo
echo "== Longline mean by age =="
grep -n "Mean by age" -A15 \
  examples/NMFS/pifsc_bigeye_tuna/level18_longline_selectivity_diagnostic/outputs/bigeye_level18_longline_prediction_decomposition.txt

echo
echo "== Longline top residuals =="
grep -n "Top residuals" -A30 \
  examples/NMFS/pifsc_bigeye_tuna/level18_longline_selectivity_diagnostic/outputs/bigeye_level18_longline_prediction_decomposition.txt
