#!/usr/bin/env bash
set -euo pipefail

echo "== FD epsilon ladder audit insertion =="
grep -n "fd_eps_ladder\\|fd_epsilon\\|cosine_analytic_fd" -A8 -B8 \
  examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/quadra/bigeye_level21_age_based_natural_mortality_diagnostic.cpp

echo
echo "== Run Level 21 with FD epsilon ladder =="
./run_bigeye_level21_age_based_m_check.sh
