#!/usr/bin/env bash
set -euo pipefail

OBJ="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/objective/bigeye_quadra_objective.hpp"

echo "== Confirm dead fleet_catch_share is gone =="
if grep -n "fleet_catch_share" "$OBJ"; then
  echo "ERROR: fleet_catch_share still present"
  exit 1
else
  echo "OK: no fleet_catch_share references remain in Level 21 objective"
fi

echo
echo "== Confirm fleet-specific Baranov terms remain =="
grep -n "catch_at_age_longline\|catch_at_age_purse_seine\|f_longline_a\|f_purse_seine_a\|const T catch_hat = is_longline" "$OBJ"

echo
echo "== Rerun Level 21 =="
./run_bigeye_level21_age_based_m_check.sh

FIT="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_fit_summary.csv"
GRAD="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_gradient_by_parameter.txt"

echo
echo "== Fit summary =="
grep -E "objective|grad_norm|converged|log_m_young_offset|log_m_old_offset" "$FIT"

echo
echo "== Top gradients =="
grep -A25 -n "Top gradients" "$GRAD"
