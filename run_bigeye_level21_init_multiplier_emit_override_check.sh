#!/usr/bin/env bash
set -euo pipefail

./run_bigeye_level21_age_based_m_check.sh

FIT="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_fit_summary.csv"
SAN="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_parameter_sanity_diagnostics.txt"

echo
echo "== Fit summary init multipliers =="
grep "init_number_multiplier_age_" "$FIT"

echo
echo "== Sanity table init multipliers =="
grep -A12 -n "Age-specific parameters" "$SAN"

echo
echo "Expected corrected fit summary:"
echo "  age1 ~1.029, age2 ~0.999, age8 ~0.724, age9 ~0.650, age10 ~0.671"
