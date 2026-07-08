#!/usr/bin/env bash
set -euo pipefail

ROOT="examples/NMFS/pifsc_bigeye_tuna/level23_longline_selectivity_smoothness_scan"

echo "== Confirm BH recruitment is synchronized in Level 23 diagnostics =="
grep -R "level23_bh_sync\|expected_recruitment\|next\[0\]" -n \
  "$ROOT/quadra/bigeye_age_structured.hpp" \
  "$ROOT/diagnostics/bigeye_longline_prediction_decomposition.hpp" \
  "$ROOT/diagnostics/bigeye_purse_seine_prediction_decomposition.hpp" \
  "$ROOT/diagnostics/bigeye_age_comp_residual_diagnostics.hpp" \
  | grep -v '\.before_' || true

echo
echo "== Confirm stale constant-R0 recruitment is gone from reports/diagnostics =="
if grep -R "next\[0\] = r0\|next\[0\] = r0 \*" -n \
  "$ROOT/quadra/bigeye_age_structured.hpp" \
  "$ROOT/diagnostics/bigeye_longline_prediction_decomposition.hpp" \
  "$ROOT/diagnostics/bigeye_purse_seine_prediction_decomposition.hpp" \
  "$ROOT/diagnostics/bigeye_age_comp_residual_diagnostics.hpp" \
  | grep -v '\.before_'; then
  echo "ERROR: stale constant-r0 recruitment remains." >&2
  exit 1
fi

echo
echo "== Run Level 23 computed-phi0 BH check =="
if [[ -x ./run_bigeye_level23_unfished_phi0_bh_check.sh ]]; then
  ./run_bigeye_level23_unfished_phi0_bh_check.sh
elif [[ -x ./run_bigeye_level23_baranov_catch_at_age_check.sh ]]; then
  ./run_bigeye_level23_baranov_catch_at_age_check.sh
else
  echo "No Level 23 run script found; compile/run manually."
fi
