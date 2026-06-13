#!/usr/bin/env bash
set -euo pipefail

FILE="core/laplace.hpp"

if [[ ! -f "$FILE" ]]; then
  echo "ERROR: $FILE not found. Run from the Quadra repo root."
  exit 1
fi

if ! git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  echo "ERROR: This does not appear to be a Git repo."
  exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
SAVE="${FILE}.saved_before_git_restore.${STAMP}"
cp "$FILE" "$SAVE"

echo "Saved current file to:"
echo "  $SAVE"
echo

echo "Restoring $FILE from Git HEAD..."
git checkout HEAD -- "$FILE"

echo
echo "Checking for bad diagnostic leftovers..."
if grep -n 'timing_hdot_end\|timing_hdot_start\|return grad;' "$FILE"; then
  echo
  echo "ERROR: Bad identifiers still exist after Git restore."
  echo "This means HEAD itself contains the bad patch."
  echo "Run:"
  echo "  git log --oneline -- core/laplace.hpp | head"
  exit 2
fi

echo "No bad diagnostic leftovers found."
echo
echo "Checking header guard closure:"
grep -n '^#ifndef QUADRA_LAPLACE_HPP\|^#define QUADRA_LAPLACE_HPP\|^#endif' "$FILE" | tail -10 || true

echo
echo "Done. Rebuild without diagnostics now."
