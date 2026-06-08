#!/usr/bin/env bash
set -euo pipefail

src="examples/state_space_surplus_production/laplace_state_space_surplus_latent_tridiagonal.cpp"
dst="examples/state_space_surplus_production/laplace_state_space_surplus_latent_runtime.cpp"

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
    "Quadra latent-state tridiagonal Laplace fixed-theta benchmark",
    "Quadra latent-state tridiagonal persistent runtime benchmark",
)

marker = "EvalResult latent_tridiagonal_laplace_eval(const ss::Data& data,"
if marker not in s:
    raise SystemExit("Could not find latent_tridiagonal_laplace_eval marker")

runtime_code = '''
class LatentTridiagonalLaplaceRuntime {
 public:
  LatentTridiagonalLaplaceRuntime(const ss::Data& data,
                                  const ss::Parameters& par)
      : data_(data),
        par_(par),
        xhat_(Eigen::VectorXd::Zero(
            static_cast<int>(data.index_observed.size() - 1))),
        initialized_(false) {}

  EvalResult evaluate_cold() {
    xhat_ = optimize_x(data_, par_);
    initialized_ = true;
    return evaluate_at_xhat();
  }

  EvalResult evaluate_warm() {
    if (!initialized_) {
      return evaluate_cold();
    }

    // Phase 2 persists xhat/backend first. Next patch will expose
    // optimize_x_from(start) so this becomes a true warm start.
    xhat_ = optimize_x(data_, par_);
    initialized_ = true;
    return evaluate_at_xhat();
  }

  EvalResult evaluate_cached_no_solve() {
    if (!initialized_) {
      throw std::runtime_error("runtime cache used before initialization");
    }
    return evaluate_at_xhat();
  }

  const char* backend_name() const {
    return backend_ ? backend_->name() : "uninitialized";
  }

  const quadra::laplace::BackendRecommendation& recommendation() const {
    return recommendation_;
  }

 private:
  EvalResult evaluate_at_xhat() {
    EvalResult out;

    out.joint = joint_x(data_, par_, xhat_);
    out.grad_norm = fd_grad_x(data_, par_, xhat_).norm();

    const Eigen::SparseMatrix<double> H =
        fd_tridiagonal_hessian_xx(data_, par_, xhat_);

    out.nnz = static_cast<int>(H.nonZeros());

    if (!backend_) {
      backend_ =
          quadra::laplace::CreateLaplaceBackendForHessian(
              H,
              &recommendation_);
    }

    backend_->factorize(H);

    if (!backend_->is_spd()) {
      throw std::runtime_error("Laplace backend reported non-SPD Hessian");
    }

    out.logdet = backend_->logdet();

    const double n_x = static_cast<double>(xhat_.size());
    out.correction = 0.5 * out.logdet -
                     0.5 * n_x * std::log(2.0 * M_PI);
    out.objective = out.joint + out.correction;

    return out;
  }

  const ss::Data& data_;
  const ss::Parameters& par_;
  Eigen::VectorXd xhat_;
  bool initialized_;

  quadra::laplace::BackendRecommendation recommendation_;
  std::unique_ptr<quadra::laplace::LaplaceBackend> backend_;
};

'''

s = s.replace(marker, runtime_code + "\n" + marker, 1)

old_loop = '''  EvalResult last = latent_tridiagonal_laplace_eval(data, par);

  const auto t0 = Clock::now();
  for (int r = 0; r < reps; ++r) {
    last = latent_tridiagonal_laplace_eval(data, par);
  }
  const auto t1 = Clock::now();

  const double total_ms = ms_between(t0, t1);
  const double avg_ms = total_ms / static_cast<double>(reps);

  std::cout << "Quadra latent-state tridiagonal Laplace fixed-theta benchmark\n";
  std::cout << "reps = " << reps << "\n";
  std::cout << "objective = " << last.objective << "\n";
  std::cout << "joint = " << last.joint << "\n";
  std::cout << "logdet = " << last.logdet << "\n";
  std::cout << "correction = " << last.correction << "\n";
  std::cout << "grad_norm = " << last.grad_norm << "\n";
  std::cout << "nnz = " << last.nnz << "\n";
  std::cout << "total_ms = " << total_ms << "\n";
  std::cout << "avg_ms = " << avg_ms << "\n";
'''

new_loop = '''  LatentTridiagonalLaplaceRuntime runtime(data, par);

  EvalResult last;

  const auto cold0 = Clock::now();
  last = runtime.evaluate_cold();
  const auto cold1 = Clock::now();
  const double cold_ms = ms_between(cold0, cold1);

  const auto warm0 = Clock::now();
  for (int r = 0; r < reps; ++r) {
    last = runtime.evaluate_warm();
  }
  const auto warm1 = Clock::now();
  const double warm_total_ms = ms_between(warm0, warm1);
  const double warm_avg_ms = warm_total_ms / static_cast<double>(reps);

  const auto cached0 = Clock::now();
  for (int r = 0; r < reps; ++r) {
    last = runtime.evaluate_cached_no_solve();
  }
  const auto cached1 = Clock::now();
  const double cached_total_ms = ms_between(cached0, cached1);
  const double cached_avg_ms = cached_total_ms / static_cast<double>(reps);

  std::cout << "Quadra latent-state tridiagonal persistent runtime benchmark\n";
  std::cout << "reps = " << reps << "\n";
  std::cout << "backend = " << runtime.backend_name() << "\n";
  std::cout << "recommendation = "
            << quadra::laplace::ToString(runtime.recommendation().backend)
            << "\n";
  std::cout << "objective = " << last.objective << "\n";
  std::cout << "joint = " << last.joint << "\n";
  std::cout << "logdet = " << last.logdet << "\n";
  std::cout << "correction = " << last.correction << "\n";
  std::cout << "grad_norm = " << last.grad_norm << "\n";
  std::cout << "nnz = " << last.nnz << "\n";
  std::cout << "cold_ms = " << cold_ms << "\n";
  std::cout << "warm_total_ms = " << warm_total_ms << "\n";
  std::cout << "warm_avg_ms = " << warm_avg_ms << "\n";
  std::cout << "cached_total_ms = " << cached_total_ms << "\n";
  std::cout << "cached_avg_ms = " << cached_avg_ms << "\n";
'''

if old_loop not in s:
    raise SystemExit("Could not find original benchmark loop")

s = s.replace(old_loop, new_loop, 1)
p.write_text(s)
PY

cat > run_state_space_surplus_latent_runtime_phase2.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-20}"
CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/examples

set -x
"${CXX}" ${CXXFLAGS}   -Iexternal/Eigen   -Iexternal/LBFGSpp/include   -I.   examples/state_space_surplus_production/laplace_state_space_surplus_latent_runtime.cpp   -o build/examples/laplace_state_space_surplus_latent_runtime

./build/examples/laplace_state_space_surplus_latent_runtime "$REPS"
EOF

chmod +x run_state_space_surplus_latent_runtime_phase2.sh

cat <<'EOF'

Installed Phase 2 state-space latent runtime example.

Run:
  ./run_state_space_surplus_latent_runtime_phase2.sh 20

Expected:
  backend = tridiagonal
  objective = -10.642176

Note:
  This persists backend and xhat.
  The warm path is not a true warm start yet because optimize_x() still starts from zero.
  Next patch will expose optimize_x_from(start).

EOF
