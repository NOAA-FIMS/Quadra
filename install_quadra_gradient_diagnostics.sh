#!/usr/bin/env bash
set -euo pipefail

# Instrument Quadra's Laplace/exact-gradient path without changing behavior.
# Run from the Quadra repository root.

TARGET="core/laplace.hpp"
if [[ ! -f "$TARGET" ]]; then
  echo "ERROR: $TARGET not found. Run this from the Quadra repo root." >&2
  exit 1
fi

BACKUP="$TARGET.bak.gradient_diagnostics.$(date +%Y%m%d_%H%M%S)"
cp "$TARGET" "$BACKUP"
echo "Backed up $TARGET -> $BACKUP"

# Add a tiny diagnostic helper include if needed.
if ! grep -q '#include <iostream>' "$TARGET"; then
  awk 'NR==1{print; print "#include <iostream>"; next} {print}' "$TARGET" > "$TARGET.tmp"
  mv "$TARGET.tmp" "$TARGET"
fi

if grep -q 'QUADRA_GRADIENT_DIAGNOSTIC' "$TARGET"; then
  echo "Diagnostic hooks already appear to be installed; leaving file unchanged."
  exit 0
fi

cat > /tmp/quadra_diag_block.txt <<'DIAG'

#if defined(QUADRA_GRADIENT_DIAGNOSTIC)
#define QUADRA_DIAG_VEC(name, v)                                                   \
  do {                                                                             \
    std::cerr << "[quadra gradient diagnostic] " << name                           \
              << " size=" << (v).size()                                           \
              << " norm=" << (v).norm()                                           \
              << " values=" << (v).transpose() << "\n";                           \
  } while (false)
#define QUADRA_DIAG_MAT(name, m)                                                   \
  do {                                                                             \
    std::cerr << "[quadra gradient diagnostic] " << name                           \
              << " rows=" << (m).rows()                                           \
              << " cols=" << (m).cols()                                           \
              << " norm=" << (m).norm() << "\n";                                  \
  } while (false)
#define QUADRA_DIAG_SCALAR(name, x)                                                \
  do {                                                                             \
    std::cerr << "[quadra gradient diagnostic] " << name << " = " << (x) << "\n"; \
  } while (false)
#else
#define QUADRA_DIAG_VEC(name, v) do {} while (false)
#define QUADRA_DIAG_MAT(name, m) do {} while (false)
#define QUADRA_DIAG_SCALAR(name, x) do {} while (false)
#endif
DIAG

# Insert macro block after includes / before first namespace-ish content.
awk '
  BEGIN { inserted=0 }
  /^#include / { print; next }
  inserted==0 { system("cat /tmp/quadra_diag_block.txt"); inserted=1 }
  { print }
' "$TARGET" > "$TARGET.tmp"
mv "$TARGET.tmp" "$TARGET"

echo "Installed diagnostic macros in $TARGET."
echo
cat <<'NEXT'
Next manual step:
  Add these calls inside the exact/profile Laplace gradient function, immediately after each value is computed:

  QUADRA_DIAG_VEC("grad_u", grad_u);
  QUADRA_DIAG_VEC("grad_theta", grad_theta);
  QUADRA_DIAG_MAT("H_u_theta", H_u_theta);
  QUADRA_DIAG_MAT("du_dtheta", du_dtheta);
  QUADRA_DIAG_VEC("implicit correction", grad_u.transpose() * du_dtheta);
  QUADRA_DIAG_VEC("logdet_grad", logdet_grad);
  QUADRA_DIAG_VEC("final analytic grad", grad);

Compile with:
  -DQUADRA_GRADIENT_DIAGNOSTIC

Then paste the diagnostic block output.
NEXT
