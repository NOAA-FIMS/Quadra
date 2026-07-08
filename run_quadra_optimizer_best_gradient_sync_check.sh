#!/usr/bin/env bash
set -euo pipefail

echo "== Optimizer best-point / fixed-gradient synchronization audit =="
grep -n "best_fx = res.value\|best_grad = grad\|best_available = true\|Eigen::VectorXd selected_grad\|selected_grad = fun\.best\|selected_grad = fun\.last_grad\|result.fixed_gradient.assign(selected_grad\|result.x = selected_x" -n core/optimizer.hpp

echo
echo "== Stale fixed-gradient assignment audit =="
if grep -n "result.fixed_gradient.assign(fun.last_grad" core/optimizer.hpp; then
  echo "ERROR: fixed_gradient still copies fun.last_grad" >&2
  exit 1
else
  echo "OK: fixed_gradient uses selected_grad, not fun.last_grad."
fi

echo
echo "== Build/run Level 21 check =="
if [[ -x ./run_bigeye_level21_age_based_m_check.sh ]]; then
  ./run_bigeye_level21_age_based_m_check.sh
else
  echo "NOTE: ./run_bigeye_level21_age_based_m_check.sh not found/executable; run your normal build manually."
fi
