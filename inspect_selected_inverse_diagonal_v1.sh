#!/usr/bin/env bash
set -euo pipefail

echo "== Selected inverse diagonal utility =="
grep -n "selected_inverse_diagonal_from_spd_hessian\\|SelectedInverseDiagonalResult\\|SimplicialLDLT" \
  core/uncertainty/selected_inverse_diagonal.hpp

echo
echo "== Test source =="
grep -n "inverse diag\\|selected_inverse_diagonal_test_passed\\|expected" \
  tests/test_selected_inverse_diagonal.cpp

echo
echo "== Run test =="
./run_selected_inverse_diagonal_test.sh
