#!/usr/bin/env bash
set -euo pipefail

TARGET="core/laplace.hpp"
if [[ ! -f "$TARGET" ]]; then
  echo "ERROR: $TARGET not found. Run this from the Quadra repo root." >&2
  exit 1
fi

latest_backup=$(ls -t core/laplace.hpp.bak.gradient_diagnostics.* 2>/dev/null | head -1 || true)
if [[ -z "$latest_backup" ]]; then
  echo "ERROR: No gradient diagnostic backup found: core/laplace.hpp.bak.gradient_diagnostics.*" >&2
  echo "Nothing changed." >&2
  exit 1
fi

save_current="core/laplace.hpp.bad_gradient_diagnostics.$(date +%Y%m%d_%H%M%S)"
cp "$TARGET" "$save_current"
cp "$latest_backup" "$TARGET"

echo "Saved current broken file -> $save_current"
echo "Restored $TARGET from -> $latest_backup"
echo
echo "Now rebuild without diagnostic changes first. If that succeeds, run:"
echo "  grep -n \"exact\|gradient\|logdet_grad\|du_dtheta\|H_u_theta\|return\" core/laplace.hpp | tail -120"
echo "and paste the relevant region so we can place diagnostics in the correct scope."
