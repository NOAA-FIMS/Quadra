#!/usr/bin/env bash
set -euo pipefail

LVL="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic"

echo "== Confirm Level 21 BH recruitment in objective =="
grep -n "unfished_spawning_biomass\|phi0\|spawning_biomass\|expected_recruitment\|next\[0\]" \
  "$LVL/objective/bigeye_quadra_objective.hpp" | head -80

echo
echo "== Stale constant recruitment audit =="
if grep -R "next\[0\] = r0\|next\[0\] = r0 \* std::exp\|next\[0\] = r0 \* exp_t" -n \
  "$LVL" --include='*.hpp' --include='*.cpp' | grep -v '\.before_'; then
  echo "WARNING: stale constant-R0 recruitment references remain above."
else
  echo "OK: no stale constant-R0 recruitment references found."
fi

echo
echo "== Build/run Level 21 age-based M check =="
if [[ -x ./run_bigeye_level21_age_based_m_check.sh ]]; then
  ./run_bigeye_level21_age_based_m_check.sh
else
  echo "No ./run_bigeye_level21_age_based_m_check.sh found; compile with your existing Level 21 runner."
fi

echo
echo "== Level 21 fit summary =="
grep -E "^(objective|grad_norm|converged|log_m_young_offset|log_m_old_offset)," \
  "$LVL/outputs/bigeye_level21_fit_summary.csv" || true
