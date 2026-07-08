#!/usr/bin/env bash
set -euo pipefail

ROOT="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic"
SRC="${ROOT}/quadra/bigeye_level21_directional_descent_check.cpp"
BIN="${ROOT}/quadra/bigeye_level21_directional_descent_check"
WF="examples/NMFS/pifsc_bigeye_tuna/workflow"

mkdir -p "${WF}"

echo "== Ensure current Level 21 free-M outputs exist =="
./run_bigeye_level21_age_based_m_check.sh

echo
echo "== Build actual directional descent checker =="
CXX="${CXX:-clang++}"
${CXX} -std=c++17 -O3 -I. -Icore -Iexternal/eigen -Iexternal -Iexternal/LBFGSpp/include \
  "${SRC}" \
  -o "${BIN}" || {
    echo
    echo "Compile failed. This probably means the objective constructor/API differs."
    echo "Inspect the normal driver construction and mirror it in:"
    echo "  ${SRC}"
    echo
    echo "Useful grep:"
    echo "  grep -n \"BigeyeQuadraObjective\\|observations\\|load\" ${ROOT}/quadra/bigeye_level21_age_based_natural_mortality_diagnostic.cpp"
    exit 1
  }

echo
echo "== Run actual directional descent checker =="
"${BIN}"

echo
echo "== Preview =="
sed -n '1,120p' "${WF}/bigeye_level21_actual_directional_descent.txt"
