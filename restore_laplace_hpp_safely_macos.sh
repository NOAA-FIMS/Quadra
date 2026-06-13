#!/usr/bin/env bash
set -euo pipefail

FILE="core/laplace.hpp"

if [[ ! -f "$FILE" ]]; then
  echo "ERROR: $FILE not found. Run this from the Quadra repo root."
  exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
BROKEN_COPY="${FILE}.broken_after_bad_cleanup.${STAMP}"
cp "$FILE" "$BROKEN_COPY"
echo "Saved current file to: $BROKEN_COPY"
echo

echo "Searching for local backups..."

BACKUPS_FILE="$(mktemp)"
find core -maxdepth 2 -type f \
  \( -name 'laplace.hpp*bak*' \
     -o -name 'laplace.hpp*backup*' \
     -o -name 'laplace.hpp*before*' \
     -o -name 'laplace.hpp*pre*' \
     -o -name 'laplace.hpp.*' \) \
  ! -name "$(basename "$BROKEN_COPY")" \
  -print 2>/dev/null | sort -r > "$BACKUPS_FILE"

if [[ -s "$BACKUPS_FILE" ]]; then
  echo "Candidate backups:"
  nl -ba "$BACKUPS_FILE"
  echo
  CHOSEN="$(head -n 1 "$BACKUPS_FILE")"
  echo "Restoring newest candidate:"
  echo "  $CHOSEN"
  cp "$CHOSEN" "$FILE"
else
  echo "No local backup found."
  echo
  if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
    echo "Falling back to Git restore of $FILE."
    echo "Your current version was saved above before restore."
    git checkout -- "$FILE"
  else
    echo "ERROR: Not in a Git repo and no local backup found."
    rm -f "$BACKUPS_FILE"
    exit 1
  fi
fi

rm -f "$BACKUPS_FILE"

echo
echo "Post-restore sanity checks:"
echo "  #ifndef count: $(grep -c '^#ifndef QUADRA_LAPLACE_HPP' "$FILE" || true)"
echo "  #define count: $(grep -c '^#define QUADRA_LAPLACE_HPP' "$FILE" || true)"
echo "  #endif count:  $(grep -c '^#endif' "$FILE" || true)"
echo
echo "Remaining bad diagnostic identifiers, if any:"
grep -n 'timing_hdot_end\|timing_hdot_start\|return grad;' "$FILE" || true
echo
echo "Done. Now rebuild without diagnostics."
