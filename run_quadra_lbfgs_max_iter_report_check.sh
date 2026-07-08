#!/usr/bin/env bash
set -euo pipefail

echo "== LBFGS max-iteration / tolerance report audit =="
grep -n "param.max_iterations\|quadra_requested_tol_met\|requested gradient tolerance\|configured max-iteration field" -A4 -B4 core/optimizer.hpp

echo
echo "== Build/run Level 21 age-based M check =="
./run_bigeye_level21_age_based_m_check.sh
