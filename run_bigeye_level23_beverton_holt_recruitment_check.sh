#!/usr/bin/env bash
set -euo pipefail

OBJ="examples/NMFS/pifsc_bigeye_tuna/level23_longline_selectivity_smoothness_scan/objective/bigeye_quadra_objective.hpp"

echo "== Confirm Level 23 Beverton-Holt recruitment formulation =="
grep -n "default_maturity_at_age\|steepness\|phi0\|spawning_biomass\|expected_recruitment\|next\[0\].*rec_dev\|r0 .* rec_dev" "${OBJ}" || true

echo
echo "== Confirm constant-R0 recruitment is gone =="
if grep -n "next\[0\] = r0 \* exp_t(rec_dev)" "${OBJ}"; then
  echo "ERROR: constant-R0 recruitment remains" >&2
  exit 1
fi

echo
echo "== Run Level 23 BH recruitment check =="
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
echo "== Level 21 vs Level 23 BH comparison =="
grep -E "objective|grad_norm|converged" \
  examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_fit_summary.csv \
  examples/NMFS/pifsc_bigeye_tuna/level23_longline_selectivity_smoothness_scan/outputs/bigeye_level23_fit_summary.csv || true

echo
echo "== Level 23 prior/residual summary =="
grep -A20 -n "Prior penalty by block" \
  examples/NMFS/pifsc_bigeye_tuna/level23_longline_selectivity_smoothness_scan/outputs/bigeye_level23_parameter_sanity_diagnostics.txt || true
grep -A8 -n "Fleet summary" \
  examples/NMFS/pifsc_bigeye_tuna/level23_longline_selectivity_smoothness_scan/outputs/bigeye_level23_age_comp_residual_diagnostics.txt || true
