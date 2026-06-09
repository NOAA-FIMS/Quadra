#!/usr/bin/env bash
set -euo pipefail
cpp="examples/opakapaka_projection/opakapaka_projection.cpp"

echo "== Reporting markers =="
grep -n "QUADRA_LEVEL1_UNCERTAINTY_REPORTING_V3\\|write_uncertainty_summary_csv\\|write_projection_uncertainty_csv" "$cpp"

echo
echo "== Run existing Opakapaka runner =="
./run_opakapaka_projection.sh

echo
echo "== New outputs =="
ls -1 examples/opakapaka_projection/outputs | grep -E 'uncertainty|covariance|correlation|standard_errors|confidence|derived|runtime' || true

echo
echo "== uncertainty_summary.csv =="
cat examples/opakapaka_projection/outputs/uncertainty_summary.csv

echo
echo "== standard_errors.csv =="
cat examples/opakapaka_projection/outputs/standard_errors.csv
