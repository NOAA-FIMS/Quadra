#!/usr/bin/env bash
set -euo pipefail

./inspect_showcase_compile_flags.sh

echo
echo "== Running O3 clean-driver Pollock showcase =="
./run_pollock_driver_showcase_report.sh

echo
echo "== Fit summary =="
cat examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_fit_summary.csv

echo
echo "== Report health excerpt =="
grep -n "Executive Summary\|Overall status\|Confidence\|Optimization quality\|Gradient quality\|Conditioning" -A8 -B2 \
  examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_analysis.md || true
