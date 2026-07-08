#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples

clang++ -std=c++17 -O3 -I"external/eigen/" \
  -DWALLEYE_POLLOCK_HUU_DIAGNOSTICS \
  -DWALLEYE_POLLOCK_HUU_MATRIX_DUMP \
  -DWALLEYE_POLLOCK_HUU_PATTERN_COMPARE \
  -DWALLEYE_POLLOCK_HUU_BAND_SUMMARY \
  -DWALLEYE_POLLOCK_HUU_BANDLIMIT_DIAGNOSTIC \
  -DWALLEYE_POLLOCK_HUU_THRESHOLD_DIAGNOSTIC \
  -DWALLEYE_POLLOCK_FIXED_HESSIAN_DIAGNOSTICS \
  -DWALLEYE_POLLOCK_LAPLACE_STRUCTURE_REPORT \
  -DWALLEYE_POLLOCK_FUNCTIONAL_ANALYSIS_REPORT \
  -DWALLEYE_POLLOCK_PARAMETER_GEOMETRY \
  -DQUADRA_LBFGS_GRAD_TOL=1.0e-2 \
  -DWALLEYE_POLLOCK_RANDOM_RECRUITMENT_COUNT=20 \
  examples/NMFS/afsc_walleye_pollock/quadra/walleye_pollock.cpp \
  examples/NMFS/afsc_walleye_pollock/quadra/walleye_pollock_adgraph_global.cpp \
  -o build/examples/pollock_functional_analysis_report

build/examples/pollock_functional_analysis_report

echo
echo "Parameter geometry section:"
awk '/Parameter Geometry/{flag=1} /Gradient Volatility/{flag=0} flag' \
  examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_functional_analysis_report.txt

echo
echo "Full functional analysis report:"
cat examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_functional_analysis_report.txt
