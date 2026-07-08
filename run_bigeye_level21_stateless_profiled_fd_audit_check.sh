#!/usr/bin/env bash
set -euo pipefail

echo "== Audit stateless profiled evaluator insertion =="
grep -n "level21_stateless_profiled_eval_v1\|stateless_profiled_eval\|const double fp =\|const double f_minus =" \
  -A12 -B8 \
  examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/quadra/bigeye_level21_age_based_natural_mortality_diagnostic.cpp

echo
echo "== Build/run Level 21 age-based M check with stateless FD audit =="
./run_bigeye_level21_age_based_m_check.sh
