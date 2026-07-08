#!/usr/bin/env bash
set -euo pipefail

ROOT="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic"

echo "== Audit remaining scalar-M diagnostic references =="
if grep -R "std::exp(-m)\|1.0 - std::exp(-m)\|const double m = level21_m_helpers::adult_m\|const double adult_m = 0.45\|const double log_m = std::log(0.45)" -n \
  "$ROOT/diagnostics" \
  --include='*.hpp' --include='*.cpp' \
  | grep -v '\.before_' | grep -v 'before_' >/tmp/level21_stale_m_refs.txt; then
  echo "WARNING: stale scalar-M diagnostic references remain:"
  cat /tmp/level21_stale_m_refs.txt
else
  echo "OK: no stale scalar-M references in Level 21 diagnostics."
fi

echo
echo "== Confirm age-specific M used in diagnostic initial equilibrium =="
grep -R "std::exp(-m_at_age\[static_cast<std::size_t>(a - 1)\]\|1.0 - std::exp(-m_at_age\[static_cast<std::size_t>(kAges - 1)\]" -n \
  "$ROOT/diagnostics" \
  --include='*.hpp' \
  | grep -v '\.before_' | grep -v 'before_' || true

echo
echo "== Build/run Level 21 age-based M check =="
./run_bigeye_level21_age_based_m_check.sh

echo
echo "== Level 21 key rows after diagnostic age-specific M sync =="
grep -E "^(objective|grad_norm|converged|log_m_young_offset|log_m_old_offset)," \
  "$ROOT/outputs/bigeye_level21_fit_summary.csv" || true

echo
echo "== Diagnostic residual summary =="
grep -n -A8 "Fleet summary" "$ROOT/outputs/bigeye_level21_age_comp_residual_diagnostics.txt" || true
