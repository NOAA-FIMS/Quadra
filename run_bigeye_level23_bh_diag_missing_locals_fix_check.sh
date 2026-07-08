#!/usr/bin/env bash
set -euo pipefail

ROOT="examples/NMFS/pifsc_bigeye_tuna/level23_longline_selectivity_smoothness_scan"

echo "== BH helper namespace openings =="
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
  echo "ERROR: expected exactly one BH helper namespace opening, found $n" >&2
  exit 1
fi

echo
echo "== Remaining BH local-symbol references =="
grep -R "const auto weight = default_weight_at_age\|const auto maturity = default_maturity_at_age\|const double phi0\|expected_recruitment\|next\[0\]" -n \
  "$ROOT/diagnostics/bigeye_longline_prediction_decomposition.hpp" \
  "$ROOT/diagnostics/bigeye_purse_seine_prediction_decomposition.hpp" \
  "$ROOT/diagnostics/bigeye_age_comp_residual_diagnostics.hpp" \
  | grep -v '\.before_' || true

echo
echo "== Run previous full checker =="
./run_bigeye_level23_bh_diag_compile_fix_check.sh
