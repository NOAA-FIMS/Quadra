#!/usr/bin/env bash
set -euo pipefail

src="examples/age_structured_recruitment/benchmark_age_structured_no_plus_flat_band.cpp"
dst="examples/age_structured_recruitment/benchmark_age_structured_no_plus_flat_band_timing.cpp"

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
    'Quadra no-plus age-structured flat-band analytic Laplace benchmark',
    'Quadra no-plus age-structured flat-band timing diagnostic'
)

old_struct = '''struct EvalResult {
  double marginal = 0.0;
  double joint = 0.0;
  double logdet = 0.0;
  double grad_norm = 0.0;
  int nnz = 0;
};'''
new_struct = '''struct EvalResult {
  double marginal = 0.0;
  double joint = 0.0;
  double logdet = 0.0;
  double grad_norm = 0.0;
  int nnz = 0;
  double optimize_ms = 0.0;
  double final_eval_ms = 0.0;
  double logdet_ms = 0.0;
  int newton_iterations = 0;
  int eval_all_calls = 0;
  int line_search_eval_calls = 0;
};

struct TimingStats {
  int newton_iterations = 0;
  int eval_all_calls = 0;
  int line_search_eval_calls = 0;
};'''
if old_struct not in s:
    raise SystemExit('missing EvalResult struct')
s = s.replace(old_struct, new_struct, 1)

s = s.replace(
    'Eigen::VectorXd optimize_x(const Data& data, const Parameters& par) {',
    'Eigen::VectorXd optimize_x(const Data& data, const Parameters& par, TimingStats* stats = nullptr) {',
    1
)
s = s.replace(
    '  double f = eval_all_flat_band(data, par, x).objective;',
    '  if (stats != nullptr) ++stats->eval_all_calls;\n  double f = eval_all_flat_band(data, par, x).objective;',
    1
)
s = s.replace(
    '    const EvalAll e = eval_all_flat_band(data, par, x);',
    '    if (stats != nullptr) {\n      ++stats->newton_iterations;\n      ++stats->eval_all_calls;\n    }\n    const EvalAll e = eval_all_flat_band(data, par, x);',
    1
)
s = s.replace(
    '      const double f_candidate = eval_all_flat_band(data, par, candidate).objective;',
    '      if (stats != nullptr) {\n        ++stats->eval_all_calls;\n        ++stats->line_search_eval_calls;\n      }\n      const double f_candidate = eval_all_flat_band(data, par, candidate).objective;',
    1
)

start = s.index('EvalResult eval_laplace(const Data& data, const Parameters& par) {')
end = s.index('\n\n}  // namespace', start)
new_eval = '''EvalResult eval_laplace(const Data& data, const Parameters& par) {
  EvalResult out;
  TimingStats stats;

  const auto opt0 = Clock::now();
  const Eigen::VectorXd xhat = optimize_x(data, par, &stats);
  const auto opt1 = Clock::now();

  const auto eval0 = Clock::now();
  const EvalAll e = eval_all_flat_band(data, par, xhat);
  const auto eval1 = Clock::now();

  const auto log0 = Clock::now();
  const double logdet = sparse_logdet_ldlt(e.hessian);
  const auto log1 = Clock::now();

  out.optimize_ms = ms_between(opt0, opt1);
  out.final_eval_ms = ms_between(eval0, eval1);
  out.logdet_ms = ms_between(log0, log1);
  out.newton_iterations = stats.newton_iterations;
  out.eval_all_calls = stats.eval_all_calls;
  out.line_search_eval_calls = stats.line_search_eval_calls;

  out.joint = e.objective;
  out.grad_norm = e.gradient.norm();
  out.nnz = static_cast<int>(e.hessian.nonZeros());
  out.logdet = logdet;

  const double n_x = static_cast<double>(xhat.size());
  out.marginal = out.joint + 0.5 * out.logdet -
                 0.5 * n_x * std::log(2.0 * M_PI);

  return out;
}'''
s = s[:start] + new_eval + s[end:]

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
            << std::setw(14) << "avg_ms"
            << std::setw(14) << "opt_ms"
            << std::setw(14) << "eval_ms"
            << std::setw(14) << "logdet_ms"
            << std::setw(10) << "newton"
            << std::setw(10) << "evals"
            << std::setw(10) << "ls_eval"
            << std::setw(10) << "nnz"
            << "\\n";'''
s = s.replace(old_header, new_header, 1)

old_row = '''    std::cout << std::setw(8) << n
              << std::setw(14) << last.marginal
              << std::setw(14) << last.joint
              << std::setw(14) << last.logdet
              << std::setw(14) << last.nnz
              << std::setw(14) << last.grad_norm
              << std::setw(14) << avg_ms
              << "\\n";'''
new_row = '''    std::cout << std::setw(8) << n
              << std::setw(14) << last.marginal
              << std::setw(14) << avg_ms
              << std::setw(14) << last.optimize_ms
              << std::setw(14) << last.final_eval_ms
              << std::setw(14) << last.logdet_ms
              << std::setw(10) << last.newton_iterations
              << std::setw(10) << last.eval_all_calls
              << std::setw(10) << last.line_search_eval_calls
              << std::setw(10) << last.nnz
              << "\\n";'''
s = s.replace(old_row, new_row, 1)

p.write_text(s)
PY

cat > run_quadra_age_structured_no_plus_flat_band_timing.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-5}"
LENGTHS="${2:-25,50,100,250,500,1000}"
AGES="${3:-10}"

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/examples

set -x
"${CXX}" ${CXXFLAGS} \
  -Iexternal/Eigen \
  examples/age_structured_recruitment/benchmark_age_structured_no_plus_flat_band_timing.cpp \
  -o build/examples/benchmark_age_structured_no_plus_flat_band_timing

./build/examples/benchmark_age_structured_no_plus_flat_band_timing "$REPS" "$LENGTHS" "$AGES"
EOF

chmod +x run_quadra_age_structured_no_plus_flat_band_timing.sh

cat <<'EOF'

Installed flat-band timing diagnostic.

Run:
  ./run_quadra_age_structured_no_plus_flat_band_timing.sh 5 25,50,100,250,500,1000 10

EOF
