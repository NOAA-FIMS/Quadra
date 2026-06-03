#!/usr/bin/env bash
set -euo pipefail

# install_age_structured_laplace_cache_demo_v1.sh
#
# Adds a small persistent-cache evaluator demo for the no-plus flat-band
# age-structured benchmark.
#
# Purpose:
#   Demonstrate the intended Quadra architecture:
#     - initialize once
#     - cache u_hat
#     - reuse u_hat as the starting point for later evaluations
#     - evaluate Laplace at cached u_hat cheaply
#
# This is a demo/scaffold, not yet wired into core LaplaceEvaluator.

src="examples/age_structured_recruitment/benchmark_age_structured_no_plus_flat_band.cpp"
dst="examples/age_structured_recruitment/laplace_cache_age_structured_no_plus_demo.cpp"

if [[ ! -f "$src" ]]; then
  echo "ERROR: missing $src"
  echo "Run install_age_structured_no_plus_flat_band_benchmark_v1.sh first."
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
    'Quadra no-plus age-structured persistent Laplace cache demo'
)

# Rename optimizer to optimize_x_from with initial x.
s = s.replace(
    'Eigen::VectorXd optimize_x(const Data& data, const Parameters& par) {\n  Eigen::VectorXd x = Eigen::VectorXd::Zero(data.n_years);',
    'Eigen::VectorXd optimize_x_from(const Data& data, const Parameters& par, const Eigen::VectorXd& initial_x) {\n  Eigen::VectorXd x = initial_x;',
    1
)

# Add default wrapper before sparse_logdet.
marker = '\n\ndouble sparse_logdet_ldlt'
wrapper = '''
Eigen::VectorXd optimize_x(const Data& data, const Parameters& par) {
  return optimize_x_from(data, par, Eigen::VectorXd::Zero(data.n_years));
}

'''
if marker not in s:
    raise SystemExit('missing sparse_logdet marker')
s = s.replace(marker, '\n\n' + wrapper + 'double sparse_logdet_ldlt', 1)

# Add eval_laplace_at_xhat and cache class before eval_laplace.
marker2 = 'EvalResult eval_laplace(const Data& data, const Parameters& par) {'
cache_code = '''
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

class PersistentLaplaceCache {
 public:
  PersistentLaplaceCache(const Data& data, const Parameters& par)
      : data_(data), par_(par), initialized_(false),
        xhat_(Eigen::VectorXd::Zero(data.n_years)) {}

  EvalResult evaluate_cold() {
    xhat_ = optimize_x(data_, par_);
    initialized_ = true;
    return eval_laplace_at_xhat(data_, par_, xhat_);
  }

  EvalResult evaluate_warm() {
    if (!initialized_) {
      return evaluate_cold();
    }
    xhat_ = optimize_x_from(data_, par_, xhat_);
    return eval_laplace_at_xhat(data_, par_, xhat_);
  }

  EvalResult evaluate_cached_no_solve() const {
    if (!initialized_) {
      throw std::runtime_error("PersistentLaplaceCache used before initialization");
    }
    return eval_laplace_at_xhat(data_, par_, xhat_);
  }

  const Eigen::VectorXd& xhat() const { return xhat_; }

 private:
  const Data& data_;
  const Parameters& par_;
  bool initialized_;
  Eigen::VectorXd xhat_;
};

'''
if marker2 not in s:
    raise SystemExit('missing eval_laplace marker')
s = s.replace(marker2, cache_code + marker2, 1)

# Replace benchmark header and loop with cache demo output.
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
            << std::setw(14) << "warm_ms"
            << std::setw(14) << "cached_ms"
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
    PersistentLaplaceCache cache(data, par);

    EvalResult last;

    const auto cold0 = Clock::now();
    last = cache.evaluate_cold();
    const auto cold1 = Clock::now();
    const double cold_ms = ms_between(cold0, cold1);

    const auto warm0 = Clock::now();
    for (int r = 0; r < reps; ++r) {
      last = cache.evaluate_warm();
    }
    const auto warm1 = Clock::now();
    const double warm_ms = ms_between(warm0, warm1) / static_cast<double>(reps);

    const auto cached0 = Clock::now();
    for (int r = 0; r < reps; ++r) {
      last = cache.evaluate_cached_no_solve();
    }
    const auto cached1 = Clock::now();
    const double cached_ms = ms_between(cached0, cached1) / static_cast<double>(reps);

    std::cout << std::setw(8) << n
              << std::setw(14) << last.marginal
              << std::setw(14) << cold_ms
              << std::setw(14) << warm_ms
              << std::setw(14) << cached_ms
              << std::setw(14) << last.grad_norm
              << std::setw(10) << last.nnz
              << "\\n";
  }
'''
if old_loop not in s:
    raise SystemExit('missing benchmark loop')
s = s.replace(old_loop, new_loop, 1)

p.write_text(s)
PY

cat > run_age_structured_laplace_cache_demo.sh <<'EOF'
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
  examples/age_structured_recruitment/laplace_cache_age_structured_no_plus_demo.cpp \
  -o build/examples/laplace_cache_age_structured_no_plus_demo

./build/examples/laplace_cache_age_structured_no_plus_demo "$REPS" "$LENGTHS" "$AGES"
EOF

chmod +x run_age_structured_laplace_cache_demo.sh

cat <<'EOF'

Installed persistent Laplace cache demo.

Run:
  ./run_age_structured_laplace_cache_demo.sh 10 25,50,100,250,500,1000 10

EOF
