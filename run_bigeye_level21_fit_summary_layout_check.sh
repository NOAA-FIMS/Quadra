#!/usr/bin/env bash
set -euo pipefail

echo "== Build/run Level 21 after fit-summary layout fix =="
./run_bigeye_level21_age_based_m_check.sh

FIT="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_fit_summary.csv"
SAN="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_parameter_sanity_diagnostics.txt"

echo
echo "== Fit summary initial-number rows =="
grep "init_number_multiplier_age_" "$FIT"

echo
echo "== Parameter sanity initial-number rows =="
grep -A15 -n "Age-specific parameters" "$SAN"

echo
echo "== Quick check: age 8/9/10 fit-summary multipliers should now be around sanity values, not 40x =="
grep -E "init_number_multiplier_age_(8|9|10)," "$FIT"
