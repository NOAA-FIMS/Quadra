#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples
mkdir -p examples/NMFS/afsc_walleye_pollock/outputs

clang++ -std=c++17 -O3 -DNDEBUG -I"external/eigen/" \
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
  -DWALLEYE_POLLOCK_MARKDOWN_REPORT \
  -DQUADRA_LBFGS_GRAD_TOL=1.0e-4 \
  -DWALLEYE_POLLOCK_RANDOM_RECRUITMENT_COUNT=20 \
  examples/NMFS/afsc_walleye_pollock/quadra/drivers/run_pollock_showcase.cpp \
  examples/NMFS/afsc_walleye_pollock/quadra/walleye_pollock_adgraph_global.cpp \
  -o build/examples/pollock_driver_showcase_report

build/examples/pollock_driver_showcase_report

echo
echo "Markdown report:"
echo "  examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_analysis.md"

echo
echo "===== POLLOCK DRIVER SHOWCASE REPORT ====="
cat examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_analysis.md
