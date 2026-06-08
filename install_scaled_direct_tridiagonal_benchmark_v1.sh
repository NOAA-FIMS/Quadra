#!/usr/bin/env bash
set -euo pipefail

src="examples/state_space_surplus_production/benchmark_latent_tridiagonal_scaled.cpp"
dst="examples/state_space_surplus_production/benchmark_latent_tridiagonal_direct_scaled.cpp"

if [[ ! -f "$src" ]]; then
  echo "ERROR: missing $src"
  exit 1
fi

cp "$src" "$dst"

python3 - "$dst" <<'PY'
import sys
from pathlib import Path

p = Path(sys.argv[1])
s = p.read_text()

s = s.replace(
    'Quadra scaled latent-state tridiagonal Laplace benchmark',
    'Quadra scaled direct latent-state tridiagonal Laplace benchmark'
)

marker = 'double sparse_logdet_ldlt(const Eigen::SparseMatrix<double>& H) {'
if marker not in s:
    raise SystemExit('Could not find sparse_logdet_ldlt marker')

helpers = (
    '\nstruct TridiagonalValues {\n'
    '  Eigen::VectorXd diag;\n'
    '  Eigen::VectorXd offdiag;\n'
    '};\n\n'
    'TridiagonalValues fd_tridiagonal_values_xx(const ss::Data& data,\n'
    '                                           const ss::Parameters& par,\n'
    '                                           const Eigen::VectorXd& x) {\n'
    '  const int n = static_cast<int>(x.size());\n'
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
    '    out.diag[j] = (gp[j] - gm[j]) / (2.0 * h);\n'
    '    if (j + 1 < n) {\n'
    '      out.offdiag[j] = (gp[j + 1] - gm[j + 1]) / (2.0 * h);\n'
    '    }\n'
    '  }\n\n'
    '  return out;\n'
    '}\n\n'
    'double logdet_tridiagonal_values_ldlt(const TridiagonalValues& H) {\n'
    '  const int n = static_cast<int>(H.diag.size());\n'
    '  if (n == 0) return 0.0;\n'
    '  double d_prev = H.diag[0];\n'
    '  if (!(d_prev > 0.0)) {\n'
    '    throw std::runtime_error("Tridiagonal value Hessian is not positive definite");\n'
    '  }\n'
    '  double logdet = std::log(d_prev);\n'
    '  for (int i = 1; i < n; ++i) {\n'
    '    const double e = H.offdiag[i - 1];\n'
    '    const double d = H.diag[i] - (e * e) / d_prev;\n'
    '    if (!(d > 0.0)) {\n'
    '      throw std::runtime_error("Tridiagonal value Hessian is not positive definite");\n'
    '    }\n'
    '    logdet += std::log(d);\n'
    '    d_prev = d;\n'
    '  }\n'
    '  return logdet;\n'
    '}\n\n'
)
s = s.replace(marker, helpers + marker, 1)

old = ('  const Eigen::SparseMatrix<double> H = fd_tridiagonal_hessian_xx(data, par, xhat);\n'
       '  out.nnz = static_cast<int>(H.nonZeros());\n'
       '  out.logdet = sparse_logdet_ldlt(H);')
new = ('  const TridiagonalValues H = fd_tridiagonal_values_xx(data, par, xhat);\n'
       '  out.nnz = static_cast<int>(H.diag.size()) +\n'
       '            2 * static_cast<int>(H.offdiag.size());\n'
       '  out.logdet = logdet_tridiagonal_values_ldlt(H);')
if old not in s:
    raise SystemExit('Could not find sparse-H eval block')
s = s.replace(old, new, 1)
p.write_text(s)
PY

cat > run_quadra_direct_vs_tmb_scaled_fixed_theta_benchmark.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-10}"
LENGTHS="${2:-25,50,100,250,500,1000}"

mkdir -p build/examples

echo "== Quadra scaled direct latent-state tridiagonal =="
set -x
c++ -std=c++17 -O2 -DNDEBUG -g \
  -Iexternal/Eigen \
  -Iexternal/LBFGSpp/include \
  -Iexamples/state_space_surplus_production \
  -Iexamples/surplus_production \
  examples/state_space_surplus_production/benchmark_latent_tridiagonal_direct_scaled.cpp \
  -o build/examples/benchmark_latent_tridiagonal_direct_scaled
./build/examples/benchmark_latent_tridiagonal_direct_scaled "$REPS" "$LENGTHS"
set +x

echo
echo "== TMB scaled AD/Laplace =="
if [[ -x ./run_tmb_scaled_state_space_surplus_benchmark.sh ]]; then
  ./run_tmb_scaled_state_space_surplus_benchmark.sh "$REPS" "$LENGTHS"
else
  ./run_quadra_vs_tmb_scaled_fixed_theta_benchmark.sh "$REPS" "$LENGTHS" | sed -n '/== TMB scaled AD\/Laplace ==/,$p'
fi
EOF

chmod +x run_quadra_direct_vs_tmb_scaled_fixed_theta_benchmark.sh

cat <<'EOF'

Installed scaled direct tridiagonal benchmark.

Run:
  ./run_quadra_direct_vs_tmb_scaled_fixed_theta_benchmark.sh 10 25,50,100,250,500,1000

EOF
