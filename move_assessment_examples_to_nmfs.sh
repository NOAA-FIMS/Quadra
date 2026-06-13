#!/usr/bin/env bash
set -euo pipefail

OLD="examples/NMFS/sefsc_red_snapper"
NEW="examples/NMFS/sefsc_red_snapper"
README="examples/NMFS/README.md"

if [[ ! -d "$OLD" && ! -d "$NEW" ]]; then
  echo "ERROR: neither $OLD nor $NEW exists. Run from repo root."
  exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
CHANGED_LIST="nmfs_example_move_changed_files_${STAMP}.txt"

echo "Preparing NMFS example reorganization..."

mkdir -p examples/NMFS

if [[ -d "$OLD" ]]; then
  if [[ -d "$NEW" ]]; then
    echo "ERROR: $NEW already exists while $OLD also exists."
    echo "Resolve manually to avoid overwriting."
    exit 1
  fi

  echo "Moving:"
  echo "  $OLD"
  echo "to:"
  echo "  $NEW"
  mv "$OLD" "$NEW"
else
  echo "$NEW already exists; skipping directory move."
fi

cat > "$README" <<'EOF'
# NMFS Assessment Examples

This directory contains fisheries stock assessment examples implemented with
Quadra.

These examples are application-oriented and are separated from smaller framework
examples so that the repository clearly distinguishes between:

- core Quadra demonstrations
- fisheries assessment model applications
- validation and comparison studies

## Current examples

### SEFSC Red Snapper

Path:

```text
examples/NMFS/sefsc_red_snapper
```

This example includes:

- age-structured population dynamics
- recruitment deviations as random effects
- Laplace approximation
- exact gradient validation
- comparison against a TMB implementation

The Red Snapper example is currently treated as a completed validation model for
Quadra's exact Laplace machinery.
EOF

echo "Updating path references..."

find . \
  -path "./.git" -prune -o \
  -type f \
  ! -name "*.o" \
  ! -name "*.so" \
  ! -name "*.dylib" \
  ! -name "*.dll" \
  ! -name "*.a" \
  ! -name "*.png" \
  ! -name "*.jpg" \
  ! -name "*.jpeg" \
  ! -name "*.pdf" \
  ! -name "*.zip" \
  ! -name "*.tar" \
  ! -name "*.gz" \
  ! -name "*.bak" \
  ! -name "*.backup" \
  ! -name "*.saved*" \
  ! -name "*.broken*" \
  -print0 |
while IFS= read -r -d '' file; do
  if grep -Iq . "$file"; then
    perl -0pi -e 's#examples/NMFS/sefsc_red_snapper#examples/NMFS/sefsc_red_snapper#g' "$file"
  fi
done

echo "Collecting changed files..."
if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  git status --short | tee "$CHANGED_LIST"
else
  find examples/NMFS -maxdepth 3 -type f | sort | tee "$CHANGED_LIST"
fi

echo
echo "Done."
echo
echo "Suggested checks:"
echo "  git status --short"
echo "  grep -R \"examples/NMFS/sefsc_red_snapper\" -n . --exclude-dir=.git"
echo
echo "Suggested Red Snapper build from repo root:"
echo '  clang++ -std=c++17 -g -I"external/eigen/" \'
echo '    examples/NMFS/sefsc_red_snapper/quadra/red_snapper_quadra_fit.cpp \'
echo '    examples/NMFS/sefsc_red_snapper/quadra/red_snapper_adgraph_global.cpp'
echo
echo "Changed-file list saved to:"
echo "  $CHANGED_LIST"
