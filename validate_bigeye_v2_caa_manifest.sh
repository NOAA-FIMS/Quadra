#!/usr/bin/env bash
set -euo pipefail

MANIFEST="examples/NMFS/pifsc_bigeye_tuna/v2/caa_manifest.yml"
BASE="examples/NMFS/pifsc_bigeye_tuna/v2/architecture/packages"

fail=0

check_package() {
  local key="$1"
  local dir="$2"

  local manifest_pkg meta_pkg
  manifest_pkg=$(awk -v key="$key" '
    $0 ~ "^  " key ":" {in_block=1; next}
    in_block && /^  [a-z_]+:/ {in_block=0}
    in_block && /package:/ {print $2; exit}
  ' "$MANIFEST")

  meta_pkg=$(awk -F': ' '/^name:/ {print $2}' "$BASE/$dir/package.meta")

  if [[ "$manifest_pkg" != "$meta_pkg" ]]; then
    echo "FAIL: $key package mismatch: manifest=$manifest_pkg meta=$meta_pkg"
    fail=1
  else
    echo "OK: $key -> $manifest_pkg"
  fi
}

check_package life_history life_history
check_package population population
check_package movement movement
check_package fleet fleet
check_package observation observation
check_package likelihood likelihood

if [[ "$fail" -ne 0 ]]; then
  exit 1
fi

echo
echo "PASSED: CAA manifest validates against package metadata"
