#!/usr/bin/env bash
set -euo pipefail

ROOT="examples/NMFS/pifsc_bigeye_tuna/level23_longline_selectivity_smoothness_scan"

echo "== Duplicate helper audit =="
grep -R "namespace level23_bh_sync {" -n \
  "$ROOT/quadra/bigeye_age_structured.hpp" \
  "$ROOT/diagnostics/bigeye_longline_prediction_decomposition.hpp" \
  "$ROOT/diagnostics/bigeye_purse_seine_prediction_decomposition.hpp" \
  "$ROOT/diagnostics/bigeye_age_comp_residual_diagnostics.hpp" \
  | grep -v '\.before_' || true

n=$(grep -R "namespace level23_bh_sync {" -n \
  "$ROOT/quadra/bigeye_age_structured.hpp" \
  "$ROOT/diagnostics/bigeye_longline_prediction_decomposition.hpp" \
  "$ROOT/diagnostics/bigeye_purse_seine_prediction_decomposition.hpp" \
  "$ROOT/diagnostics/bigeye_age_comp_residual_diagnostics.hpp" \
  | grep -v '\.before_' | wc -l | tr -d ' ')

if [[ "$n" != "1" ]]; then
  echo "ERROR: expected exactly one BH helper namespace after cleanup, found $n" >&2
  exit 1
fi

echo
echo "== Required local symbols =="
grep -R "const auto weight = default_weight_at_age\|const auto maturity = default_maturity_at_age\|const double phi0\|expected_recruitment\|next\[0\]" -n \
  "$ROOT/diagnostics/bigeye_longline_prediction_decomposition.hpp" \
  "$ROOT/diagnostics/bigeye_purse_seine_prediction_decomposition.hpp" \
  "$ROOT/diagnostics/bigeye_age_comp_residual_diagnostics.hpp" \
  | grep -v '\.before_' || true

echo
echo "== Stale constant recruitment audit =="
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
echo "== Build/run Level 23 computed-phi0 BH check =="
if [[ -x ./run_bigeye_level23_unfished_phi0_bh_check.sh ]]; then
  ./run_bigeye_level23_unfished_phi0_bh_check.sh
else
  echo "No run_bigeye_level23_unfished_phi0_bh_check.sh found."
fi
