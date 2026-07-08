#!/usr/bin/env bash
set -euo pipefail

echo "== O3 build Bigeye Level 16 objective-consistent residuals =="
mkdir -p build/examples
c++ -std=c++17 -O3 -flto \
  -I. -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level16_purse_seine_age_selectivity/quadra/bigeye_level16_purse_seine_age_selectivity.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level16_purse_seine_age_selectivity/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level16_purse_seine_age_selectivity_check

echo
echo "== Run Bigeye Level 16 objective-consistent residuals =="
./build/examples/pifsc_bigeye_level16_purse_seine_age_selectivity_check

echo
echo "== Age-comp residual summary =="
grep -n "Fleet summary" -A8 \
  examples/NMFS/pifsc_bigeye_tuna/level16_purse_seine_age_selectivity/outputs/bigeye_level16_age_comp_residual_diagnostics.txt

echo
echo "== Worst residuals =="
grep -n "Compact worst age-comp residuals" -A25 \
  examples/NMFS/pifsc_bigeye_tuna/level16_purse_seine_age_selectivity/outputs/bigeye_level16_age_comp_residual_diagnostics.txt

echo
echo "== Prediction decomposition mean by age =="
grep -n "Mean by age" -A15 \
  examples/NMFS/pifsc_bigeye_tuna/level16_purse_seine_age_selectivity/outputs/bigeye_level16_purse_seine_prediction_decomposition.txt || true
