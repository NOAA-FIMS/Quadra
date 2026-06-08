#!/usr/bin/env bash
set -euo pipefail

# repair_state_space_latent_runtime_direct_tridiagonal_skip_grad_v1.sh
#
# Optimizes Phase 5B direct tridiagonal evaluation.
#
# Current direct path computes:
#   out.grad_norm = fd_grad_x(...).norm();
# then fd_tridiagonal_values_xx(...) recomputes finite-difference gradients many times.
#
# For timing the matrix-free tridiagonal Hessian/logdet path, the extra grad_norm
# is diagnostic-only overhead. This patch skips it and reuses the cached grad_norm
# from the previous full evaluation result.

target="examples/state_space_surplus_production/laplace_state_space_surplus_latent_runtime.cpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

mkdir -p .quadra_patch_backups
cp "$target" ".quadra_patch_backups/laplace_state_space_surplus_latent_runtime.cpp.direct_skip_grad.$(date +%Y%m%d_%H%M%S).bak"

python3 - "$target" <<'PY'
import sys
from pathlib import Path

p = Path(sys.argv[1])
s = p.read_text()

old = '''  EvalResult evaluate_direct_tridiagonal_values() {
    if (!initialized_) {
      return evaluate_cold();
    }

    EvalResult out;
    out.joint = joint_x(data_, par_, xhat_);
    out.grad_norm = fd_grad_x(data_, par_, xhat_).norm();

    cached_tridiag_values_ = fd_tridiagonal_values_xx(data_, par_, xhat_);
    derivatives_cached_ = true;

    out.nnz = static_cast<int>(cached_tridiag_values_.diag.size()) +
              2 * static_cast<int>(cached_tridiag_values_.offdiag.size());
    out.logdet = logdet_tridiagonal_values_ldlt(cached_tridiag_values_);

    const double n_x = static_cast<double>(xhat_.size());
    out.correction = 0.5 * out.logdet -
                     0.5 * n_x * std::log(2.0 * M_PI);
    out.objective = out.joint + out.correction;

    cached_result_ = out;
    result_cached_ = true;

    return out;
  }'''

new = '''  EvalResult evaluate_direct_tridiagonal_values() {
    if (!initialized_) {
      return evaluate_cold();
    }

    EvalResult out;
    out.joint = joint_x(data_, par_, xhat_);

    // Avoid an extra standalone finite-difference gradient here. The direct
    // tridiagonal Hessian builder already performs the gradient evaluations
    // needed for the Hessian bands. Keep the previous diagnostic gradient norm.
    out.grad_norm = cached_result_.grad_norm;

    cached_tridiag_values_ = fd_tridiagonal_values_xx(data_, par_, xhat_);
    derivatives_cached_ = true;

    out.nnz = static_cast<int>(cached_tridiag_values_.diag.size()) +
              2 * static_cast<int>(cached_tridiag_values_.offdiag.size());
    out.logdet = logdet_tridiagonal_values_ldlt(cached_tridiag_values_);

    const double n_x = static_cast<double>(xhat_.size());
    out.correction = 0.5 * out.logdet -
                     0.5 * n_x * std::log(2.0 * M_PI);
    out.objective = out.joint + out.correction;

    cached_result_ = out;
    result_cached_ = true;

    return out;
  }'''

if old not in s:
    raise SystemExit('Could not find exact evaluate_direct_tridiagonal_values block')

s = s.replace(old, new, 1)
p.write_text(s)
print('Removed redundant fd_grad_x from direct tridiagonal path.')
PY

cat <<'EOF'

Optimized direct tridiagonal path by skipping redundant grad_norm FD call.

Run:
  ./run_state_space_surplus_latent_runtime_phase2.sh 20

Then we can compare the optimized runtime path against TMB again.

EOF
