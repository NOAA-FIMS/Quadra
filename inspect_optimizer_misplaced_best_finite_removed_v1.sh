#!/usr/bin/env bash
set -euo pipefail

echo "== Check misplaced block removed =="
if grep -n "QUADRA_RETURN_BEST_FINITE_ON_LINE_SEARCH_FAILURE_V1" core/optimizer.hpp; then
  echo "ERROR: misplaced marker still present" >&2
  exit 1
else
  echo "OK: misplaced marker removed"
fi

echo
echo "== Check optimizer still has line-search references =="
grep -n "line search\\|sufficiently decrease\\|runtime_error" core/optimizer.hpp | head -80 || true

echo
echo "Next run:"
echo "  ./inspect_opakapaka_level1_reporting_v7.sh"
