#!/usr/bin/env bash
set -euo pipefail

./inspect_bigeye_level10_remove_stale_geometry.sh

echo
echo "== O3 build Bigeye Level 10 q anchor after geometry cleanup =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level10_q_anchor/quadra/bigeye_level10_q_anchor.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level10_q_anchor/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level10_q_anchor_check

echo
echo "== Run Bigeye Level 10 q anchor =="
./build/examples/pifsc_bigeye_level10_q_anchor_check

echo
echo "== Level 10 fit summary =="
cat examples/NMFS/pifsc_bigeye_tuna/level10_q_anchor/outputs/bigeye_level10_fit_summary.csv

echo
echo "== Level 10 components =="
cat examples/NMFS/pifsc_bigeye_tuna/level10_q_anchor/outputs/bigeye_level10_objective_components.csv

echo
echo "== Level 10 recruitment diagnostics =="
sed -n '1,150p' examples/NMFS/pifsc_bigeye_tuna/level10_q_anchor/outputs/bigeye_level10_recruitment_diagnostics.txt
