#!/usr/bin/env bash
set -euo pipefail

LVL="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic"

echo "== Level 21 diagnostic phi0 context =="
grep -R "const double phi0\|beverton_holt_recruitment\|next\[0\]" -n \
  "$LVL/diagnostics" --include='*.hpp' | grep -v '\.before_' || true

echo
echo "== Build/run Level 21 BH check =="
./run_bigeye_level21_age_based_m_check.sh
