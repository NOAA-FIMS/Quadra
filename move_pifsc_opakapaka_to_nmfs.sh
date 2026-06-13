#!/usr/bin/env bash
set -euo pipefail

OLD="examples/NMFS/pifsc_opakapaka"
NEW="examples/NMFS/pifsc_opakapaka"
README="examples/NMFS/README.md"

if [[ ! -d "$OLD" && ! -d "$NEW" ]]; then
  echo "ERROR: neither $OLD nor $NEW exists. Run from repo root."
  exit 1
fi

STAMP="$(date +%Y%m%d_%H%M%S)"
CHANGED_LIST="pifsc_opakapaka_nmfs_move_changed_files_${STAMP}.txt"

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

echo "Updating active path references..."

# Update active repo files, but avoid:
# - .git
# - build artifacts
# - historical patch backups
# - previously generated changed-file inventories
# - compiled Opakapaka executable
find . \
  -path "./.git" -prune -o \
  -path "./.quadra_patch_backups" -prune -o \
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
  ! -name "nmfs_example_move_changed_files_*.txt" \
  ! -name "pifsc_opakapaka_nmfs_move_changed_files_*.txt" \
  ! -path "./examples/NMFS/pifsc_opakapaka/quadra/opakapaka_projection" \
  -print0 |
while IFS= read -r -d '' file; do
  if grep -Iq . "$file"; then
    perl -0pi -e 's#examples/NMFS/pifsc_opakapaka#examples/NMFS/pifsc_opakapaka#g' "$file"
  fi
done

echo "Updating examples/NMFS/README.md..."

if [[ ! -f "$README" ]]; then
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
EOF
fi

if ! grep -q "PIFSC Opakapaka" "$README"; then
  cat >> "$README" <<'EOF'

### PIFSC Opakapaka

Path:

```text
examples/NMFS/pifsc_opakapaka
```

This example includes:

- Pacific Islands assessment-style projection workflow
- synthetic data input
- uncertainty reporting
- derived quantities
- projection uncertainty outputs
- comparison against a TMB implementation
EOF
fi

echo "Collecting changed files..."
if git rev-parse --is-inside-work-tree >/dev/null 2>&1; then
  git status --short | tee "$CHANGED_LIST"
else
  find examples/NMFS -maxdepth 4 -type f | sort | tee "$CHANGED_LIST"
fi

echo
echo "Remaining active references to old path, excluding backups and git:"
grep -R "examples/NMFS/pifsc_opakapaka" -n . \
  --exclude-dir=.git \
  --exclude-dir=.quadra_patch_backups \
  --exclude="*.bak" \
  --exclude="*.backup" \
  --exclude="nmfs_example_move_changed_files_*.txt" \
  --exclude="pifsc_opakapaka_nmfs_move_changed_files_*.txt" || true

echo
echo "Done."
echo
echo "Suggested build/check:"
echo '  clang++ -std=c++17 -g -I"external/eigen/" \'
echo '    examples/NMFS/pifsc_opakapaka/quadra/opakapaka_projection.cpp \'
echo '    -o build/examples/NMFS/pifsc_opakapaka'
echo
echo "Changed-file list saved to:"
echo "  $CHANGED_LIST"
