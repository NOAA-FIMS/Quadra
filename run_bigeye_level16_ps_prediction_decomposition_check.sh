#!/usr/bin/env bash
set -euo pipefail

echo "== O3 build Bigeye Level 16 PS prediction decomposition =="
mkdir -p build/examples
c++ -std=c++17 -O3 -flto \
  -I. -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level16_purse_seine_age_selectivity/quadra/bigeye_level16_purse_seine_age_selectivity.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level16_purse_seine_age_selectivity/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level16_purse_seine_age_selectivity_check

echo
echo "== Run Bigeye Level 16 PS prediction decomposition =="
./build/examples/pifsc_bigeye_level16_purse_seine_age_selectivity_check

echo
echo "== Prediction decomposition preview =="
sed -n '1,140p' \
  examples/NMFS/pifsc_bigeye_tuna/level16_purse_seine_age_selectivity/outputs/bigeye_level16_purse_seine_prediction_decomposition.txt

echo
echo "== Mean by age =="
grep -n "Mean by age" -A15 \
  examples/NMFS/pifsc_bigeye_tuna/level16_purse_seine_age_selectivity/outputs/bigeye_level16_purse_seine_prediction_decomposition.txt

echo
echo "== Top residuals =="
grep -n "Top residuals" -A25 \
  examples/NMFS/pifsc_bigeye_tuna/level16_purse_seine_age_selectivity/outputs/bigeye_level16_purse_seine_prediction_decomposition.txt
