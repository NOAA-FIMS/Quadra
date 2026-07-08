#!/usr/bin/env bash
set -euo pipefail

ROOT="examples/NMFS/pifsc_bigeye_tuna"

echo "== Check for remaining fleet-average F patterns =="
grep -R "0\.5.*sel_longline.*sel_purse_seine\|T(0\.5).*sel_longline.*sel_purse_seine" -n \
  "$ROOT"/level20_longline_selectivity_regularization_scan \
  "$ROOT"/level21_age_based_natural_mortality_diagnostic \
  "$ROOT"/level23_longline_selectivity_smoothness_scan \
  --include='*.hpp' --include='*.cpp' || true

echo
echo "== Check total-selectivity F patterns =="
grep -R "total_sel.*sel_longline.*sel_purse_seine\|total_sel_prev.*sel_longline.*sel_purse_seine\|total_sel_last.*sel_longline.*sel_purse_seine" -n \
  "$ROOT"/level20_longline_selectivity_regularization_scan \
  "$ROOT"/level21_age_based_natural_mortality_diagnostic \
  "$ROOT"/level23_longline_selectivity_smoothness_scan \
  --include='*.hpp' --include='*.cpp' | head -160

echo
echo "== Rerun Level 21 age-M check =="
./run_bigeye_level21_age_based_m_check.sh
