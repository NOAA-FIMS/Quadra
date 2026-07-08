#!/usr/bin/env bash
set -euo pipefail

ROOT="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic"

echo "== Audit undeclared/scalar m remnants in Level 21 diagnostics =="
if grep -R "std::exp(-m)\|1.0 - std::exp(-m)\|z_prev = m +\|z_last = m +\|<< m <<\|m_fixed" -n \
  "$ROOT/diagnostics" \
  --include='*.hpp' --include='*.cpp' \
  | grep -v '\.before_' | grep -v 'before_' >/tmp/level21_remaining_m_refs.txt; then
  echo "WARNING: possible scalar-M remnants remain:"
  cat /tmp/level21_remaining_m_refs.txt
else
  echo "OK: no obvious scalar-M diagnostic remnants."
fi

echo
echo "== Confirm age-specific M dynamics in patched diagnostics =="
grep -R "z_prev = m_at_age\[prev\]\|z_last = m_at_age\[last\]\|m_young:\|summary,m_young" -n \
  "$ROOT/diagnostics" \
  --include='*.hpp' \
  | grep -v '\.before_' | grep -v 'before_' || true

echo
echo "== Build/run Level 21 age-based M check =="
./run_bigeye_level21_age_based_m_check.sh
