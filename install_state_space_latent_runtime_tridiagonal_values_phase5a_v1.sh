#!/usr/bin/env bash
set -euo pipefail

target="examples/state_space_surplus_production/laplace_state_space_surplus_latent_runtime.cpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

mkdir -p .quadra_patch_backups
cp "$target" ".quadra_patch_backups/laplace_state_space_surplus_latent_runtime.cpp.phase5a_tridiag_values.$(date +%Y%m%d_%H%M%S).bak"

python3 - "$target" <<'PY'
import sys
from pathlib import Path

p = Path(sys.argv[1])
s = p.read_text()

if 'struct TridiagonalValues' in s:
    print('Phase 5A tridiagonal values patch already appears installed.')
    raise SystemExit(0)

marker = 'double sparse_logdet_ldlt(const Eigen::SparseMatrix<double>& H) {'
if marker not in s:
    raise SystemExit('Could not find sparse_logdet_ldlt marker')

helpers = (
    '\nstruct TridiagonalValues {\n'
    '  Eigen::VectorXd diag;\n'
    '  Eigen::VectorXd offdiag;  // offdiag[i - 1] = H(i, i - 1)\n'
    '};\n\n'
    'TridiagonalValues fd_tridiagonal_values_xx(const ss::Data& data,\n'
    '                                           const ss::Parameters& par,\n'
    '                                           const Eigen::VectorXd& x) {\n'
    '  const int n = static_cast<int>(x.size());\n\n'
    '  TridiagonalValues out;\n'
    '  out.diag = Eigen::VectorXd::Zero(n);\n'
    '  out.offdiag = Eigen::VectorXd::Zero(std::max(0, n - 1));\n\n'
    '  for (int j = 0; j < n; ++j) {\n'
    '    const double h = 1e-5 * (1.0 + std::abs(x[j]));\n'
    '    Eigen::VectorXd xp = x;\n'
    '    Eigen::VectorXd xm = x;\n'
    '    xp[j] += h;\n'
    '    xm[j] -= h;\n\n'
    '    const Eigen::VectorXd gp = fd_grad_x(data, par, xp);\n'
    '    const Eigen::VectorXd gm = fd_grad_x(data, par, xm);\n\n'
    '    out.diag[j] = (gp[j] - gm[j]) / (2.0 * h);\n\n'
    '    if (j + 1 < n) {\n'
    '      out.offdiag[j] = (gp[j + 1] - gm[j + 1]) / (2.0 * h);\n'
    '    }\n'
    '  }\n\n'
    '  return out;\n'
    '}\n\n'
    'double logdet_tridiagonal_values_ldlt(const TridiagonalValues& H) {\n'
    '  const int n = static_cast<int>(H.diag.size());\n'
    '  if (n == 0) return 0.0;\n\n'
    '  double d_prev = H.diag[0];\n'
    '  if (!(d_prev > 0.0)) {\n'
    '    throw std::runtime_error("Tridiagonal value Hessian is not positive definite");\n'
    '  }\n\n'
    '  double logdet = std::log(d_prev);\n\n'
    '  for (int i = 1; i < n; ++i) {\n'
    '    const double e = H.offdiag[i - 1];\n'
    '    const double d = H.diag[i] - (e * e) / d_prev;\n\n'
    '    if (!(d > 0.0)) {\n'
    '      throw std::runtime_error("Tridiagonal value Hessian is not positive definite");\n'
    '    }\n\n'
    '    logdet += std::log(d);\n'
    '    d_prev = d;\n'
    '  }\n\n'
    '  return logdet;\n'
    '}\n\n'
)
s = s.replace(marker, helpers + marker, 1)

anchor = '  Eigen::SparseMatrix<double> cached_H_;\n'
if anchor not in s:
    raise SystemExit('Could not find cached_H_ member anchor; install derivative cache first')
s = s.replace(anchor, anchor + '  TridiagonalValues cached_tridiag_values_;\n', 1)

anchor = '    cached_H_ = H;\n    derivatives_cached_ = true;\n'
if anchor not in s:
    raise SystemExit('Could not find cached_H_ assignment anchor')
s = s.replace(anchor, anchor + '    cached_tridiag_values_ = fd_tridiagonal_values_xx(data_, par_, xhat_);\n', 1)

anchor = '  const char* backend_name() const {'
if anchor not in s:
    raise SystemExit('Could not find backend_name anchor')
method = (
    '  EvalResult evaluate_cached_tridiagonal_values() {\n'
    '    if (!initialized_ || !derivatives_cached_ || !result_cached_) {\n'
    '      throw std::runtime_error("runtime tridiagonal value cache used before initialization");\n'
    '    }\n\n'
    '    EvalResult out = cached_result_;\n'
    '    out.logdet = logdet_tridiagonal_values_ldlt(cached_tridiag_values_);\n\n'
    '    const double n_x = static_cast<double>(xhat_.size());\n'
    '    out.correction = 0.5 * out.logdet -\n'
    '                     0.5 * n_x * std::log(2.0 * M_PI);\n'
    '    out.objective = out.joint + out.correction;\n\n'
    '    return out;\n'
    '  }\n\n'
)
s = s.replace(anchor, method + anchor, 1)

anchor = '  const auto result_cached0 = Clock::now();\n'
if anchor not in s:
    raise SystemExit('Could not find result_cached0 timing anchor')
timing = (
    '  const auto tridiagonal_values0 = Clock::now();\n'
    '  for (int r = 0; r < reps; ++r) {\n'
    '    last = runtime.evaluate_cached_tridiagonal_values();\n'
    '  }\n'
    '  const auto tridiagonal_values1 = Clock::now();\n'
    '  const double tridiagonal_values_total_ms =\n'
    '      ms_between(tridiagonal_values0, tridiagonal_values1);\n'
    '  const double tridiagonal_values_avg_ms =\n'
    '      tridiagonal_values_total_ms / static_cast<double>(reps);\n\n'
)
s = s.replace(anchor, timing + anchor, 1)

anchor = '  std::cout << "result_cached_total_ms = " << result_cached_total_ms << "\\n";\n'
if anchor not in s:
    raise SystemExit('Could not find result_cached_total_ms output anchor')
output = (
    '  std::cout << "tridiagonal_values_total_ms = " << tridiagonal_values_total_ms << "\\n";\n'
    '  std::cout << "tridiagonal_values_avg_ms = " << tridiagonal_values_avg_ms << "\\n";\n'
)
s = s.replace(anchor, output + anchor, 1)

p.write_text(s)
print('Installed Phase 5A tridiagonal values cache.')
PY

cat <<'EOF'

Installed Phase 5A tridiagonal values cache benchmark.

Run:
  ./run_state_space_surplus_latent_runtime_phase2.sh 20

Expected additional output:
  tridiagonal_values_total_ms = ...
  tridiagonal_values_avg_ms = ...

EOF
