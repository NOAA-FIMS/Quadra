#!/usr/bin/env bash
set -euo pipefail

src="examples/age_structured_recruitment/benchmark_age_structured_no_plus_flat_band.cpp"
dst="examples/age_structured_recruitment/benchmark_age_structured_no_plus_flat_band_newton_diagnostics.cpp"

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

s = s.replace('Quadra no-plus age-structured flat-band analytic Laplace benchmark',
              'Quadra no-plus age-structured flat-band Newton diagnostic')

marker = 'Eigen::VectorXd optimize_x(const Data& data, const Parameters& par) {'
diag_struct = '''
struct SolveDiagnostics {
  int newton_iterations = 0;
  int eval_all_calls = 0;
  int line_search_eval_calls = 0;
  double final_grad_norm = 0.0;
};

'''
if marker not in s:
    raise SystemExit('missing optimize_x marker')
s = s.replace(marker, diag_struct + marker, 1)

s = s.replace('Eigen::VectorXd optimize_x(const Data& data, const Parameters& par) {',
              'Eigen::VectorXd optimize_x_diagnostic(const Data& data, const Parameters& par, SolveDiagnostics* diag = nullptr) {', 1)

s = s.replace('  double f = eval_all_flat_band(data, par, x).objective;',
              '  if (diag != nullptr) ++diag->eval_all_calls;\n  double f = eval_all_flat_band(data, par, x).objective;', 1)

s = s.replace('    const EvalAll e = eval_all_flat_band(data, par, x);',
              '    if (diag != nullptr) {\n      ++diag->newton_iterations;\n      ++diag->eval_all_calls;\n    }\n    const EvalAll e = eval_all_flat_band(data, par, x);\n    if (diag != nullptr) diag->final_grad_norm = e.gradient.norm();', 1)

s = s.replace('      const double f_candidate = eval_all_flat_band(data, par, candidate).objective;',
              '      if (diag != nullptr) {\n        ++diag->eval_all_calls;\n        ++diag->line_search_eval_calls;\n      }\n      const double f_candidate = eval_all_flat_band(data, par, candidate).objective;', 1)

marker2 = '\n\ndouble sparse_logdet_ldlt'
wrapper = '''
Eigen::VectorXd optimize_x(const Data& data, const Parameters& par) {
  return optimize_x_diagnostic(data, par, nullptr);
}

'''
if marker2 not in s:
    raise SystemExit('missing sparse_logdet marker')
s = s.replace(marker2, '\n\n' + wrapper + 'double sparse_logdet_ldlt', 1)

marker3 = 'EvalResult eval_laplace(const Data& data, const Parameters& par) {'
helper = '''
EvalResult eval_laplace_at_xhat(const Data& data,
                                const Parameters& par,
                                const Eigen::VectorXd& xhat) {
  const EvalAll e = eval_all_flat_band(data, par, xhat);

  EvalResult out;
  out.joint = e.objective;
  out.grad_norm = e.gradient.norm();
  out.nnz = static_cast<int>(e.hessian.nonZeros());
  out.logdet = sparse_logdet_ldlt(e.hessian);

  const double n_x = static_cast<double>(xhat.size());
  out.marginal = out.joint + 0.5 * out.logdet -
                 0.5 * n_x * std::log(2.0 * M_PI);

  return out;
}

'''
if marker3 not in s:
    raise SystemExit('missing eval_laplace marker')
s = s.replace(marker3, helper + marker3, 1)

# Replace header by locating known sequence compactly.
old_header = '''  std::cout << std::setw(8) << "n"
            << std::setw(14) << "objective"
            << std::setw(14) << "joint"
            << std::setw(14) << "logdet"
            << std::setw(14) << "nnz"
            << std::setw(14) << "grad_norm"
            << std::setw(14) << "avg_ms"
            << "\\n";'''
new_header = '''  std::cout << std::setw(8) << "n"
            << std::setw(14) << "objective"
            << std::setw(14) << "cold_ms"
            << std::setw(10) << "newton"
            << std::setw(10) << "evals"
            << std::setw(10) << "ls_eval"
            << std::setw(14) << "grad_norm"
            << std::setw(10) << "nnz"
            << "\\n";'''
if old_header not in s:
    raise SystemExit('missing header')
s = s.replace(old_header, new_header, 1)

old_loop = '''  for (const int n : lengths) {
    const Data data = make_data(n, n_ages, par);
    EvalResult last = eval_laplace(data, par);

    const auto t0 = Clock::now();
    for (int r = 0; r < reps; ++r) {
      last = eval_laplace(data, par);
    }
    const auto t1 = Clock::now();

    const double avg_ms = ms_between(t0, t1) / static_cast<double>(reps);

    std::cout << std::setw(8) << n
              << std::setw(14) << last.marginal
              << std::setw(14) << last.joint
              << std::setw(14) << last.logdet
              << std::setw(14) << last.nnz
              << std::setw(14) << last.grad_norm
              << std::setw(14) << avg_ms
              << "\\n";
  }
'''
new_loop = '''  for (const int n : lengths) {
    const Data data = make_data(n, n_ages, par);

    EvalResult last;
    SolveDiagnostics diag;

    const auto t0 = Clock::now();
    for (int r = 0; r < reps; ++r) {
      diag = SolveDiagnostics{};
      Eigen::VectorXd xhat = optimize_x_diagnostic(data, par, &diag);
      last = eval_laplace_at_xhat(data, par, xhat);
    }
    const auto t1 = Clock::now();

    const double cold_ms = ms_between(t0, t1) / static_cast<double>(reps);

    std::cout << std::setw(8) << n
              << std::setw(14) << last.marginal
              << std::setw(14) << cold_ms
              << std::setw(10) << diag.newton_iterations
              << std::setw(10) << diag.eval_all_calls
              << std::setw(10) << diag.line_search_eval_calls
              << std::setw(14) << last.grad_norm
              << std::setw(10) << last.nnz
              << "\\n";
  }
'''
if old_loop not in s:
    raise SystemExit('missing loop')
s = s.replace(old_loop, new_loop, 1)

p.write_text(s)
PY

cat > run_quadra_age_structured_no_plus_flat_band_newton_diagnostics.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-3}"
LENGTHS="${2:-25,50,100,250,500,1000}"
AGES="${3:-10}"

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/examples

set -x
"${CXX}" ${CXXFLAGS} \
  -Iexternal/Eigen \
  examples/age_structured_recruitment/benchmark_age_structured_no_plus_flat_band_newton_diagnostics.cpp \
  -o build/examples/benchmark_age_structured_no_plus_flat_band_newton_diagnostics

./build/examples/benchmark_age_structured_no_plus_flat_band_newton_diagnostics "$REPS" "$LENGTHS" "$AGES"
EOF

chmod +x run_quadra_age_structured_no_plus_flat_band_newton_diagnostics.sh

cat <<'EOF'

Installed flat-band Newton diagnostic.

Run:
  ./run_quadra_age_structured_no_plus_flat_band_newton_diagnostics.sh 3 25,50,100,250,500,1000 10

EOF
