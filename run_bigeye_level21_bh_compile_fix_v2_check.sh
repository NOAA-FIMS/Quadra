#!/usr/bin/env bash
set -euo pipefail

LVL="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic"

echo "== Objective BH/phi0 context =="
grep -n "maturity\|steepness\|unfished_spawning_biomass\|phi0\|spawning_biomass\|expected_recruitment\|next\[0\]" \
  "$LVL/objective/bigeye_quadra_objective.hpp" | head -140

echo
echo "== Stale constant recruitment audit =="
if grep -R "next\[0\] = r0\|next\[0\] = r0 \* std::exp\|next\[0\] = r0 \* exp_t" -n \
  "$LVL" --include='*.hpp' --include='*.cpp' | grep -v '\.before_'; then
  echo "ERROR: stale constant recruitment remains."
  exit 1
else
  echo "OK: no stale constant-R0 recruitment references found."
fi

echo
echo "== Build/run Level 21 BH check =="
./run_bigeye_level21_age_based_m_check.sh
