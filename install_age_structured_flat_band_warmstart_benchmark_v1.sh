#!/usr/bin/env bash
set -euo pipefail

src="examples/age_structured_recruitment/benchmark_age_structured_no_plus_flat_band.cpp"
dst="examples/age_structured_recruitment/benchmark_age_structured_no_plus_flat_band_warmstart.cpp"

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
    'Quadra no-plus age-structured flat-band warm-start benchmark'
)

s = s.replace(
    'Eigen::VectorXd optimize_x(const Data& data, const Parameters& par) {\n  Eigen::VectorXd x = Eigen::VectorXd::Zero(data.n_years);',
    'Eigen::VectorXd optimize_x_from(const Data& data, const Parameters& par, const Eigen::VectorXd& initial_x) {\n  Eigen::VectorXd x = initial_x;',
    1
)

marker = '\n\ndouble sparse_logdet_ldlt'
wrapper = '''
Eigen::VectorXd optimize_x(const Data& data, const Parameters& par) {
  return optimize_x_from(data, par, Eigen::VectorXd::Zero(data.n_years));
}

'''
if marker not in s:
    raise SystemExit('missing sparse_logdet marker')
s = s.replace(marker, '\n\n' + wrapper + 'double sparse_logdet_ldlt', 1)

marker2 = 'EvalResult eval_laplace(const Data& data, const Parameters& par) {'
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
if marker2 not in s:
    raise SystemExit('missing eval_laplace marker')
s = s.replace(marker2, helper + marker2, 1)

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
            << std::setw(14) << "cached_ms"
            << std::setw(14) << "warm_ms"
            << std::setw(14) << "grad_norm"
            << std::setw(10) << "nnz"
            << "\\n";'''
if old_header not in s:
    raise SystemExit('missing table header')
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

    const auto cold0 = Clock::now();
    Eigen::VectorXd xhat = optimize_x(data, par);
    EvalResult last = eval_laplace_at_xhat(data, par, xhat);
    const auto cold1 = Clock::now();
    const double cold_ms = ms_between(cold0, cold1);

    const auto cached0 = Clock::now();
    for (int r = 0; r < reps; ++r) {
      last = eval_laplace_at_xhat(data, par, xhat);
    }
    const auto cached1 = Clock::now();
    const double cached_ms = ms_between(cached0, cached1) / static_cast<double>(reps);

    const auto warm0 = Clock::now();
    for (int r = 0; r < reps; ++r) {
      xhat = optimize_x_from(data, par, xhat);
      last = eval_laplace_at_xhat(data, par, xhat);
    }
    const auto warm1 = Clock::now();
    const double warm_ms = ms_between(warm0, warm1) / static_cast<double>(reps);

    std::cout << std::setw(8) << n
              << std::setw(14) << last.marginal
              << std::setw(14) << cold_ms
              << std::setw(14) << cached_ms
              << std::setw(14) << warm_ms
              << std::setw(14) << last.grad_norm
              << std::setw(10) << last.nnz
              << "\\n";
  }
'''
if old_loop not in s:
    raise SystemExit('missing main loop')
s = s.replace(old_loop, new_loop, 1)

p.write_text(s)
PY

cat > run_quadra_age_structured_no_plus_flat_band_warmstart.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-10}"
LENGTHS="${2:-25,50,100,250,500,1000}"
AGES="${3:-10}"

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/examples

set -x
"${CXX}" ${CXXFLAGS} \
  -Iexternal/Eigen \
  examples/age_structured_recruitment/benchmark_age_structured_no_plus_flat_band_warmstart.cpp \
  -o build/examples/benchmark_age_structured_no_plus_flat_band_warmstart

./build/examples/benchmark_age_structured_no_plus_flat_band_warmstart "$REPS" "$LENGTHS" "$AGES"
EOF

chmod +x run_quadra_age_structured_no_plus_flat_band_warmstart.sh

cat <<'EOF'

Installed flat-band warm-start/cached xhat benchmark.

Run:
  ./run_quadra_age_structured_no_plus_flat_band_warmstart.sh 10 25,50,100,250,500,1000 10

Interpretation:
  cold_ms   = solve xhat from zero + final Laplace eval
  cached_ms = final Laplace eval at already-solved xhat
  warm_ms   = optimize starting from previous xhat + final Laplace eval

EOF
