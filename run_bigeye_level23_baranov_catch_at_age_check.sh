#!/usr/bin/env bash
set -euo pipefail

OBJ="examples/NMFS/pifsc_bigeye_tuna/level23_longline_selectivity_smoothness_scan/objective/bigeye_quadra_objective.hpp"

echo "== Confirm Level 23 fleet-specific Baranov catch-at-age formulation =="
grep -n "catch_at_age_longline\|catch_at_age_purse_seine\|f_longline_a\|f_purse_seine_a\|f_total_a\|fleet_catch_at_age\|catch_hat" "${OBJ}" || true

echo
echo "== Confirm old formulation is gone =="
if grep -n "total_catch_hat\|fleet_catch_share\|pred_age_comp\[i\] = n\[i\].*sel" "${OBJ}"; then
  echo "ERROR: old catch/share or vulnerable-number composition remnants remain" >&2
  exit 1
fi

echo
echo "== Run Level 23 check =="
if [[ -x ./run_bigeye_level23_longline_selectivity_smoothness_check.sh ]]; then
  ./run_bigeye_level23_longline_selectivity_smoothness_check.sh
elif [[ -x ./run_bigeye_level23_longline_selectivity_smoothness_scan.sh ]]; then
  ./run_bigeye_level23_longline_selectivity_smoothness_scan.sh
else
  echo "No known Level 23 runner found in repo root."
  find examples/NMFS/pifsc_bigeye_tuna/level23_longline_selectivity_smoothness_scan -maxdepth 3 -type f -name "run*.sh" -print | sort
  exit 1
fi

echo
echo "== Level 21 vs Level 23 aligned formulation comparison =="
grep -E "objective|grad_norm|converged" \
  examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_fit_summary.csv \
  examples/NMFS/pifsc_bigeye_tuna/level23_longline_selectivity_smoothness_scan/outputs/bigeye_level23_fit_summary.csv || true

echo
echo "== Re-run checklist if available =="
if [[ -x ./audit_bigeye_multifan_consistency_checklist.sh ]]; then
  ./audit_bigeye_multifan_consistency_checklist.sh
fi
