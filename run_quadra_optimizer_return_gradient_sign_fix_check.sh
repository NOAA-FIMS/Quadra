#!/usr/bin/env bash
set -euo pipefail

echo "== Confirm optimizer boundary gradient sign fix =="
grep -n "quadra_optimizer_return_gradient_sign_fix_v1\|grad = -to_eigen(res.grad_x)" -A12 -B4 core/optimizer.hpp

echo
echo "== Build/run Level 21 age-based M check =="
./run_bigeye_level21_age_based_m_check.sh
