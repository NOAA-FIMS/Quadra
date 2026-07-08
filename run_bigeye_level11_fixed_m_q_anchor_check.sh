#!/usr/bin/env bash
set -euo pipefail

./inspect_bigeye_level11_fixed_m_q_anchor.sh

echo
echo "== O3 build Bigeye Level 11 fixed M + q anchor =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level11_fixed_m_q_anchor/quadra/bigeye_level11_fixed_m_q_anchor.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level11_fixed_m_q_anchor/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level11_fixed_m_q_anchor_check

echo
echo "== Run Bigeye Level 11 fixed M + q anchor =="
./build/examples/pifsc_bigeye_level11_fixed_m_q_anchor_check

echo
echo "== Level 11 fit summary =="
cat examples/NMFS/pifsc_bigeye_tuna/level11_fixed_m_q_anchor/outputs/bigeye_level11_fit_summary.csv

echo
echo "== Level 11 components =="
cat examples/NMFS/pifsc_bigeye_tuna/level11_fixed_m_q_anchor/outputs/bigeye_level11_objective_components.csv

echo
echo "== Level 11 recruitment diagnostics =="
sed -n '1,150p' examples/NMFS/pifsc_bigeye_tuna/level11_fixed_m_q_anchor/outputs/bigeye_level11_recruitment_diagnostics.txt
