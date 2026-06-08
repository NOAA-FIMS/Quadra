#!/usr/bin/env bash
set -euo pipefail

# install_state_space_latent_runtime_derivative_cache_v2.sh
# Self-contained robust derivative-cache patch.

target="examples/state_space_surplus_production/laplace_state_space_surplus_latent_runtime.cpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

mkdir -p .quadra_patch_backups
cp "$target" ".quadra_patch_backups/laplace_state_space_surplus_latent_runtime.cpp.derivative_cache_v2.$(date +%Y%m%d_%H%M%S).bak"

python3 - "$target" <<'PY'
import sys
from pathlib import Path

p = Path(sys.argv[1])
s = p.read_text()

# Add derivative cache members.
if 'Eigen::SparseMatrix<double> cached_H_' not in s:
    anchor = '  bool initialized_;\n'
    if anchor not in s:
        raise SystemExit('Could not find initialized_ member anchor')
    s = s.replace(anchor, anchor + '  bool derivatives_cached_ = false;\n  Eigen::SparseMatrix<double> cached_H_;\n', 1)

# If result cache members are absent, add them too.
if 'EvalResult cached_result_' not in s:
    anchor = '  bool initialized_;\n'
    if anchor not in s:
        raise SystemExit('Could not find initialized_ member anchor for result cache')
    s = s.replace(anchor, anchor + '  bool result_cached_ = false;\n  EvalResult cached_result_;\n', 1)

# Cache H in evaluate_at_xhat after construction.
if 'cached_H_ = H;' not in s:
    old = ('    const Eigen::SparseMatrix<double> H =\n'
           '        fd_tridiagonal_hessian_xx(data_, par_, xhat_);\n\n'
           '    out.nnz = static_cast<int>(H.nonZeros());')
    new = ('    const Eigen::SparseMatrix<double> H =\n'
           '        fd_tridiagonal_hessian_xx(data_, par_, xhat_);\n\n'
           '    cached_H_ = H;\n'
           '    derivatives_cached_ = true;\n\n'
           '    out.nnz = static_cast<int>(H.nonZeros());')
    if old not in s:
        raise SystemExit('Could not find H construction block')
    s = s.replace(old, new, 1)

# Cache result in evaluate_at_xhat.
if 'cached_result_ = out;' not in s:
    old = '    out.objective = out.joint + out.correction;\n\n    return out;\n  }'
    new = ('    out.objective = out.joint + out.correction;\n\n'
           '    cached_result_ = out;\n'
           '    result_cached_ = true;\n\n'
           '    return out;\n  }')
    if old not in s:
        raise SystemExit('Could not find evaluate_at_xhat return block')
    s = s.replace(old, new, 1)

# Add evaluate_cached_result if missing.
if 'EvalResult evaluate_cached_result() const' not in s:
    anchor = '  const char* backend_name() const {'
    method = (
        '  EvalResult evaluate_cached_result() const {\n'
        '    if (!initialized_ || !result_cached_) {\n'
        '      throw std::runtime_error("runtime result cache used before initialization");\n'
        '    }\n'
        '    return cached_result_;\n'
        '  }\n\n'
    )
    if anchor not in s:
        raise SystemExit('Could not find backend_name anchor for result cache method')
    s = s.replace(anchor, method + anchor, 1)

# Add evaluate_cached_derivatives if missing.
if 'EvalResult evaluate_cached_derivatives()' not in s:
    anchor = '  const char* backend_name() const {'
    method = (
        '  EvalResult evaluate_cached_derivatives() {\n'
        '    if (!initialized_ || !derivatives_cached_ || !result_cached_) {\n'
        '      throw std::runtime_error("runtime derivative cache used before initialization");\n'
        '    }\n\n'
        '    EvalResult out = cached_result_;\n\n'
        '    if (!backend_) {\n'
        '      backend_ =\n'
        '          quadra::laplace::CreateLaplaceBackendForHessian(\n'
        '              cached_H_,\n'
        '              &recommendation_);\n'
        '    }\n\n'
        '    backend_->factorize(cached_H_);\n\n'
        '    if (!backend_->is_spd()) {\n'
        '      throw std::runtime_error("Laplace backend reported non-SPD Hessian");\n'
        '    }\n\n'
        '    out.logdet = backend_->logdet();\n'
        '    const double n_x = static_cast<double>(xhat_.size());\n'
        '    out.correction = 0.5 * out.logdet -\n'
        '                     0.5 * n_x * std::log(2.0 * M_PI);\n'
        '    out.objective = out.joint + out.correction;\n\n'
        '    cached_result_ = out;\n'
        '    result_cached_ = true;\n\n'
        '    return out;\n'
        '  }\n\n'
    )
    if anchor not in s:
        raise SystemExit('Could not find backend_name anchor for derivative method')
    s = s.replace(anchor, method + anchor, 1)

# Add derivative timing before result-cache timing.
if 'derivative_cached_total_ms' not in s:
    anchor = '  const auto result_cached0 = Clock::now();\n'
    insert = (
        '  const auto derivative_cached0 = Clock::now();\n'
        '  for (int r = 0; r < reps; ++r) {\n'
        '    last = runtime.evaluate_cached_derivatives();\n'
        '  }\n'
        '  const auto derivative_cached1 = Clock::now();\n'
        '  const double derivative_cached_total_ms =\n'
        '      ms_between(derivative_cached0, derivative_cached1);\n'
        '  const double derivative_cached_avg_ms =\n'
        '      derivative_cached_total_ms / static_cast<double>(reps);\n\n'
    )
    if anchor not in s:
        raise SystemExit('Could not find result_cached0 timing anchor')
    s = s.replace(anchor, insert + anchor, 1)

# Add derivative output lines.
if 'derivative_cached_avg_ms = ' not in s:
    anchor = '  std::cout << "result_cached_total_ms = " << result_cached_total_ms << "\\n";\n'
    insert = (
        '  std::cout << "derivative_cached_total_ms = " << derivative_cached_total_ms << "\\n";\n'
        '  std::cout << "derivative_cached_avg_ms = " << derivative_cached_avg_ms << "\\n";\n'
    )
    if anchor not in s:
        raise SystemExit('Could not find result_cached_total_ms output anchor')
    s = s.replace(anchor, insert + anchor, 1)

p.write_text(s)
print('Installed derivative cache v2.')
PY

cat <<'EOF'

Installed derivative cache v2.

Run:
  ./run_state_space_surplus_latent_runtime_phase2.sh 20

Expected additional output:
  derivative_cached_total_ms = ...
  derivative_cached_avg_ms = ...

EOF
