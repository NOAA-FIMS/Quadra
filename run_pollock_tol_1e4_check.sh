#!/usr/bin/env bash
set -euo pipefail

./inspect_pollock_tolerance_setting.sh

echo
echo "== Running Pollock polished report with tolerance 1e-4 =="
./run_pollock_polished_report.sh

echo
echo "== Post-run fit summary =="
cat examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_fit_summary.csv

echo
echo "== Post-run model health excerpt =="
grep -n "Executive Summary\|Model Health Assessment\|Optimization\|Gradient quality\|Overall status\|Confidence" -A10 -B2 \
  examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_analysis.md || true
