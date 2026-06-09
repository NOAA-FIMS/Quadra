#!/usr/bin/env bash
set -euo pipefail

echo "== Optimizer best-finite line-search return markers =="
grep -n "QUADRA_RETURN_BEST_FINITE_ON_LINE_SEARCH_FAILURE_V1\\|best finite iterate after LBFGS line-search failure\\|sufficiently decrease" core/optimizer.hpp

echo
echo "== Relevant best_finite symbols =="
grep -n "best_finite_seen\\|best_finite_x\\|best_finite_fx\\|best_finite_grad_norm\\|best_finite_u_hat\\|eigen_to_std_vector\\|fixed_grad_tol" core/optimizer.hpp | head -80

echo
echo "Now rerun:"
echo "  ./inspect_opakapaka_level1_reporting_v6.sh"
echo "or:"
echo "  ./inspect_opakapaka_level1_reporting_v7.sh"
