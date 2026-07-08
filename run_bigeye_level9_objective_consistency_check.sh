#!/usr/bin/env bash
set -euo pipefail

./inspect_bigeye_level9_objective_consistency.sh

echo
echo "== O3 build Bigeye Level 9 objective consistency =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level9_estimated_natural_mortality/quadra/bigeye_level9_estimated_natural_mortality.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level9_estimated_natural_mortality/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level9_objective_consistency_check

echo
echo "== Run Bigeye Level 9 objective consistency =="
./build/examples/pifsc_bigeye_level9_objective_consistency_check

echo
echo "== Objective consistency check =="
cat examples/NMFS/pifsc_bigeye_tuna/level9_estimated_natural_mortality/outputs/bigeye_level9_objective_consistency_check.csv

echo
echo "== Objective components =="
cat examples/NMFS/pifsc_bigeye_tuna/level9_estimated_natural_mortality/outputs/bigeye_level9_objective_components.csv || true

echo
echo "== Fit summary =="
cat examples/NMFS/pifsc_bigeye_tuna/level9_estimated_natural_mortality/outputs/bigeye_level9_fit_summary.csv
