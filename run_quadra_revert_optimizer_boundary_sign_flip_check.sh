#!/usr/bin/env bash
set -euo pipefail

echo "== Confirm optimizer boundary sign convention =="
grep -n "quadra_optimizer_boundary_gradient_convention_v2\|grad = to_eigen(res.grad_x)\|grad = -to_eigen(res.grad_x)" -A12 -B4 core/optimizer.hpp

if grep -q "grad = -to_eigen(res.grad_x)" core/optimizer.hpp; then
  echo "ERROR: boundary sign flip is still present." >&2
  exit 1
fi

echo
echo "== Build/run Level 21 age-based M check =="
./run_bigeye_level21_age_based_m_check.sh
