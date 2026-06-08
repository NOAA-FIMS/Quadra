#!/usr/bin/env bash
set -euo pipefail

# install_flat_band_newton_alpha_trace_v1.sh
# Adds per-Newton-iteration alpha/backtracking trace for the no-plus flat-band age benchmark.

src="examples/age_structured_recruitment/benchmark_age_structured_no_plus_flat_band.cpp"
dst="examples/age_structured_recruitment/benchmark_age_structured_no_plus_flat_band_alpha_trace.cpp"

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
              'Quadra no-plus age-structured flat-band Newton alpha trace')

# Add trace flag global inside namespace, after namespace opening.
s = s.replace('namespace {\n\n', 'namespace {\n\nbool g_trace_newton_alpha = false;\n\n', 1)

# Inject trace variables after gradient norm.
s = s.replace(
    '    const double gnorm = e.gradient.norm();\n    if (gnorm < 1e-8) break;',
    '    const double gnorm = e.gradient.norm();\n    if (gnorm < 1e-8) break;\n\n    const double current_f = f;',
    1
)

# Add step norm after dx solve.
s = s.replace(
    '    const Eigen::VectorXd dx = ldlt.solve(e.gradient);\n    if (ldlt.info() != Eigen::Success || !dx.allFinite()) {',
    '    const Eigen::VectorXd dx = ldlt.solve(e.gradient);\n    const double step_norm = dx.norm();\n    if (ldlt.info() != Eigen::Success || !dx.allFinite()) {',
    1
)

# Add backtrack counter.
s = s.replace(
    '    double step = 1.0;\n    bool accepted = false;\n\n    for (int ls = 0; ls < 30; ++ls) {',
    '    double step = 1.0;\n    bool accepted = false;\n    int accepted_ls = -1;\n\n    for (int ls = 0; ls < 30; ++ls) {',
    1
)

# Capture accepted ls.
s = s.replace(
    '        accepted = true;\n        break;',
    '        accepted = true;\n        accepted_ls = ls;\n        break;',
    1
)

# Print per iteration after line-search block before failure check.
s = s.replace(
    '    if (!accepted) {\n      if (gnorm < 1e-5) break;',
    '    if (g_trace_newton_alpha) {\n      std::cout << "iter " << std::setw(3) << iter\n                << " f " << std::setw(14) << current_f\n                << " gnorm " << std::setw(14) << gnorm\n                << " step_norm " << std::setw(14) << step_norm\n                << " alpha " << std::setw(14) << step\n                << " backtracks " << std::setw(4) << accepted_ls\n                << " accepted " << (accepted ? "yes" : "no")\n                << "\\n";\n    }\n\n    if (!accepted) {\n      if (gnorm < 1e-5) break;',
    1
)

# Replace main loop with trace-only selected sizes.
old_header = '''  std::cout << std::setw(8) << "n"
            << std::setw(14) << "objective"
            << std::setw(14) << "joint"
            << std::setw(14) << "logdet"
            << std::setw(14) << "nnz"
            << std::setw(14) << "grad_norm"
            << std::setw(14) << "avg_ms"
            << "\\n";'''
new_header = '''  std::cout << "Tracing Newton alpha/backtracking only.\\n";
'''
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
    std::cout << "\\n============================================================\\n";
    std::cout << "n = " << n << ", ages = " << n_ages << "\\n";
    std::cout << "============================================================\\n";

    g_trace_newton_alpha = true;
    const auto t0 = Clock::now();
    const Eigen::VectorXd xhat = optimize_x(data, par);
    const auto t1 = Clock::now();
    g_trace_newton_alpha = false;

    const EvalResult last = eval_laplace(data, par);
    std::cout << "summary objective=" << last.marginal
              << " grad_norm=" << last.grad_norm
              << " nnz=" << last.nnz
              << " solve_ms=" << ms_between(t0, t1)
              << "\\n";
  }
'''
if old_loop not in s:
    raise SystemExit('missing loop')
s = s.replace(old_loop, new_loop, 1)

p.write_text(s)
PY

cat > run_quadra_age_structured_no_plus_flat_band_alpha_trace.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

LENGTHS="${1:-500,1000}"
AGES="${2:-10}"

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/examples

set -x
"${CXX}" ${CXXFLAGS} \
  -Iexternal/Eigen \
  examples/age_structured_recruitment/benchmark_age_structured_no_plus_flat_band_alpha_trace.cpp \
  -o build/examples/benchmark_age_structured_no_plus_flat_band_alpha_trace

./build/examples/benchmark_age_structured_no_plus_flat_band_alpha_trace 1 "$LENGTHS" "$AGES"
EOF

chmod +x run_quadra_age_structured_no_plus_flat_band_alpha_trace.sh

cat <<'EOF'

Installed flat-band Newton alpha trace.

Run:
  ./run_quadra_age_structured_no_plus_flat_band_alpha_trace.sh 500,1000 10

EOF
