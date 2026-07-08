#!/usr/bin/env bash
set -euo pipefail

echo "== Static check: Level 21 fleet-specific Baranov catch-at-age =="
OBJ="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/objective/bigeye_quadra_objective.hpp"
grep -n "longline_catch_hat\|purse_seine_catch_hat\|longline_catch_at_age\|purse_seine_catch_at_age\|fleet_catch_at_age_sum\|fleet_catch_share\|total_catch_hat\|selected_numbers_sum" "$OBJ" || true

echo
echo "== Build/run Level 21 fleet-specific Baranov objective =="
./run_bigeye_level21_age_based_m_check.sh

FIT="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_fit_summary.csv"
GRAD="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_gradient_by_parameter.txt"
SAN="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_parameter_sanity_diagnostics.txt"
RES="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_age_comp_residual_diagnostics.txt"

echo
echo "== Fit summary =="
grep -E "objective|grad_norm|converged|log_m_young_offset|log_m_old_offset|init_number_multiplier_age_" "$FIT"

echo
echo "== Top gradients =="
grep -A30 -n "Top gradients" "$GRAD"

echo
echo "== Parameter sanity prior blocks =="
grep -A18 -n "Prior penalty by block" "$SAN"

echo
echo "== Residual fleet summary =="
grep -A8 -n "Fleet summary" "$RES"
