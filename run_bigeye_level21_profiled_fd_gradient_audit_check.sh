#!/usr/bin/env bash
set -euo pipefail

echo "== Audit inserted profiled FD gradient block =="
grep -n "Profiled FD Gradient Audit\|cosine_analytic_fd\|analytic_minus_fd\|fd_eps" -A8 -B8 \
  examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/quadra/bigeye_level21_age_based_natural_mortality_diagnostic.cpp

echo
echo "== Build/run Level 21 age-based M check with profiled FD gradient audit =="
./run_bigeye_level21_age_based_m_check.sh
