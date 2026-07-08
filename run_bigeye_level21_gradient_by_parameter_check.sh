#!/usr/bin/env bash
set -euo pipefail

./run_bigeye_level21_age_based_m_check.sh

echo
echo "== Level 21 top fixed-effect gradients =="
grep -n "Summary\|max_abs_fixed_gradient\|Top gradients" -A35 \
  examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_gradient_by_parameter.txt || true

echo
echo "== Compact top 15 fixed gradients =="
awk -F, 'NR==1 || (NR>1 && NR<=16)' \
  examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/outputs/bigeye_level21_gradient_by_parameter.csv || true
