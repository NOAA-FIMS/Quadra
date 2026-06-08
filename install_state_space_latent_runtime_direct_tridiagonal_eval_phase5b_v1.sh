#!/usr/bin/env bash
set -euo pipefail

# install_state_space_latent_runtime_direct_tridiagonal_eval_phase5b_v1.sh
#
# Phase 5B:
#   Add a direct tridiagonal-values evaluation path for cached xhat.
#
# New method:
#   evaluate_direct_tridiagonal_values()
#
# It computes:
#   joint
#   grad_norm
#   TridiagonalValues diag/offdiag
#   matrix-free tridiagonal LDLT logdet
#
# without constructing SparseMatrix H.

target="examples/state_space_surplus_production/laplace_state_space_surplus_latent_runtime.cpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

mkdir -p .quadra_patch_backups
cp "$target" ".quadra_patch_backups/laplace_state_space_surplus_latent_runtime.cpp.phase5b_direct_tridiag.$(date +%Y%m%d_%H%M%S).bak"

python3 - "$target" <<'PY'
import sys
from pathlib import Path

p = Path(sys.argv[1])
s = p.read_text()

if 'evaluate_direct_tridiagonal_values()' in s:
    print('Phase 5B direct tridiagonal path already appears installed.')
    raise SystemExit(0)

if 'struct TridiagonalValues' not in s:
    raise SystemExit('Missing TridiagonalValues helpers; run Phase 5A installer first')

anchor = '  const char* backend_name() const {'
if anchor not in s:
    raise SystemExit('Could not find backend_name anchor')

method = (
    '  EvalResult evaluate_direct_tridiagonal_values() {\n'
    '    if (!initialized_) {\n'
    '      return evaluate_cold();\n'
    '    }\n\n'
    '    EvalResult out;\n'
    '    out.joint = joint_x(data_, par_, xhat_);\n'
    '    out.grad_norm = fd_grad_x(data_, par_, xhat_).norm();\n\n'
    '    cached_tridiag_values_ = fd_tridiagonal_values_xx(data_, par_, xhat_);\n'
    '    derivatives_cached_ = true;\n\n'
    '    out.nnz = static_cast<int>(cached_tridiag_values_.diag.size()) +\n'
    '              2 * static_cast<int>(cached_tridiag_values_.offdiag.size());\n'
    '    out.logdet = logdet_tridiagonal_values_ldlt(cached_tridiag_values_);\n\n'
    '    const double n_x = static_cast<double>(xhat_.size());\n'
    '    out.correction = 0.5 * out.logdet -\n'
    '                     0.5 * n_x * std::log(2.0 * M_PI);\n'
    '    out.objective = out.joint + out.correction;\n\n'
    '    cached_result_ = out;\n'
    '    result_cached_ = true;\n\n'
    '    return out;\n'
    '  }\n\n'
)
s = s.replace(anchor, method + anchor, 1)

# Add timing after cached_no_solve timing and before derivative cached timing.
anchor = '  const auto derivative_cached0 = Clock::now();\n'
if anchor not in s:
    raise SystemExit('Could not find derivative_cached0 timing anchor; run derivative cache v2 first')

timing = (
    '  const auto direct_tridiagonal0 = Clock::now();\n'
    '  for (int r = 0; r < reps; ++r) {\n'
    '    last = runtime.evaluate_direct_tridiagonal_values();\n'
    '  }\n'
    '  const auto direct_tridiagonal1 = Clock::now();\n'
    '  const double direct_tridiagonal_total_ms =\n'
    '      ms_between(direct_tridiagonal0, direct_tridiagonal1);\n'
    '  const double direct_tridiagonal_avg_ms =\n'
    '      direct_tridiagonal_total_ms / static_cast<double>(reps);\n\n'
)
s = s.replace(anchor, timing + anchor, 1)

# Add output before derivative cache output.
anchor = '  std::cout << "derivative_cached_total_ms = " << derivative_cached_total_ms << "\\n";\n'
if anchor not in s:
    raise SystemExit('Could not find derivative_cached output anchor')

output = (
    '  std::cout << "direct_tridiagonal_total_ms = " << direct_tridiagonal_total_ms << "\\n";\n'
    '  std::cout << "direct_tridiagonal_avg_ms = " << direct_tridiagonal_avg_ms << "\\n";\n'
)
s = s.replace(anchor, output + anchor, 1)

p.write_text(s)
print('Installed Phase 5B direct tridiagonal evaluation path.')
PY

cat <<'EOF'

Installed Phase 5B direct tridiagonal evaluation path.

Run:
  ./run_state_space_surplus_latent_runtime_phase2.sh 20

Expected additional output:
  direct_tridiagonal_total_ms = ...
  direct_tridiagonal_avg_ms = ...

Interpretation:
  cached_avg_ms              = rebuild sparse H + backend logdet
  direct_tridiagonal_avg_ms  = compute diag/offdiag + matrix-free logdet
  derivative_cached_avg_ms   = reuse sparse H + backend logdet
  tridiagonal_values_avg_ms  = reuse diag/offdiag + matrix-free logdet

EOF
