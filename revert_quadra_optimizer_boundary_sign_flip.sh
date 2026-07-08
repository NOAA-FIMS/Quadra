#!/usr/bin/env bash
set -euo pipefail

TARGET="core/optimizer.hpp"
STAMP="$(date +%Y%m%d_%H%M%S)"
BACKUP="${TARGET}.before_revert_optimizer_boundary_sign_flip.${STAMP}"

if [[ ! -f "$TARGET" ]]; then
  echo "ERROR: target not found: $TARGET" >&2
  exit 1
fi

cp "$TARGET" "$BACKUP"
echo "backup: $BACKUP"

python3 - <<'PY'
from pathlib import Path

path = Path("core/optimizer.hpp")
s = path.read_text()

old = """    // quadra_optimizer_return_gradient_sign_fix_v1
    //
    // The Level 21 profiled-FD audit showed the analytic fixed-effect vector
    // returned by laplace_eval_at_u_star_persistent_structured is opposite the
    // finite-difference gradient of the profiled objective. LBFGSpp expects the
    // true objective gradient, so normalize the sign at the optimizer boundary.
    //
    // This keeps the lower-level Laplace pieces untouched while making:
    //   x - alpha * grad
    // a decreasing direction for the profiled objective.
    grad = -to_eigen(res.grad_x);
"""

new = """    // quadra_optimizer_boundary_gradient_convention_v2
    //
    // Do not flip this sign at the optimizer boundary.
    //
    // Empirical Level 21 check:
    //   - flipping this to -to_eigen(res.grad_x) caused LBFGS to stall at the
    //     initial objective (~4255) with |grad| ~1004.
    //   - the unflipped convention descends to the Level 21 basin
    //     (~1731.585) before line-search precision limits.
    //
    // The remaining sign mismatch is now treated as an audit/convention issue
    // to isolate lower in the Laplace-gradient stack, not as an optimizer-boundary
    // fix.
    grad = to_eigen(res.grad_x);
"""

if old in s:
    s = s.replace(old, new, 1)
elif "grad = -to_eigen(res.grad_x);" in s:
    s = s.replace("    grad = -to_eigen(res.grad_x);\n", new, 1)
else:
    print("No boundary sign flip found; leaving optimizer gradient assignment unchanged.")

path.write_text(s)
print(f"patched: {path}")
PY

cat > run_quadra_revert_optimizer_boundary_sign_flip_check.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

echo "== Confirm optimizer boundary sign convention =="
grep -n "quadra_optimizer_boundary_gradient_convention_v2\|grad = to_eigen(res.grad_x)\|grad = -to_eigen(res.grad_x)" -A12 -B4 core/optimizer.hpp

if grep -q "grad = -to_eigen(res.grad_x)" core/optimizer.hpp; then
  echo "ERROR: boundary sign flip is still present." >&2
  exit 1
fi

echo
echo "== Build/run Level 21 age-based M check =="
./run_bigeye_level21_age_based_m_check.sh
SH

chmod +x run_quadra_revert_optimizer_boundary_sign_flip_check.sh

echo
echo "Created:"
echo "  ./run_quadra_revert_optimizer_boundary_sign_flip_check.sh"
echo
echo "Run:"
echo "  ./run_quadra_revert_optimizer_boundary_sign_flip_check.sh"
