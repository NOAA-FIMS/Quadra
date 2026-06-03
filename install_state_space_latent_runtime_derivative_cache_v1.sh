#!/usr/bin/env bash
set -euo pipefail

# install_state_space_latent_runtime_derivative_cache_v1.sh
#
# Phase 4:
#   Add an explicit derivative artifact cache to the runtime example.
#
# This is different from result_cached:
#   result_cached_result() returns the full last EvalResult directly.
#
# New derivative-cache mode:
#   evaluate_cached_derivatives()
#
# It reuses:
#   - cached joint
#   - cached grad_norm
#   - cached H
#   - cached nnz
#
# but still exercises:
#   - backend->factorize(H)
#   - backend->logdet()
#   - correction/objective assembly
#
# This isolates backend/factorization cost from FD derivative reconstruction.
#
# Production design note:
#   invalidate_derivatives() should be called whenever theta or xhat changes.

target="examples/state_space_surplus_production/laplace_state_space_surplus_latent_runtime.cpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

mkdir -p .quadra_patch_backups
cp "$target" ".quadra_patch_backups/laplace_state_space_surplus_latent_runtime.cpp.derivative_cache.$(date +%Y%m%d_%H%M%S).bak"

python3 - "$target" <<'PY'
import sys
from pathlib import Path

p = Path(sys.argv[1])
s = p.read_text()

if "evaluate_cached_derivatives()" in s:
    print("Derivative cache patch already appears installed.")
    raise SystemExit(0)

# Add derivative-cache members.
old_members = '''  bool initialized_;
  bool result_cached_ = false;
  EvalResult cached_result_;

  quadra::laplace::BackendRecommendation recommendation_;'''

new_members = '''  bool initialized_;
  bool result_cached_ = false;
  bool derivatives_cached_ = false;
  EvalResult cached_result_;
  Eigen::SparseMatrix<double> cached_H_;

  quadra::laplace::BackendRecommendation recommendation_;'''

if old_members not in s:
    raise SystemExit("Could not find member block with result cache")

s = s.replace(old_members, new_members, 1)

# In evaluate_at_xhat, cache H when it is built.
old_h_block = '''    const Eigen::SparseMatrix<double> H =
        fd_tridiagonal_hessian_xx(data_, par_, xhat_);

    out.nnz = static_cast<int>(H.nonZeros());'''

new_h_block = '''    const Eigen::SparseMatrix<double> H =
        fd_tridiagonal_hessian_xx(data_, par_, xhat_);

    cached_H_ = H;
    derivatives_cached_ = true;

    out.nnz = static_cast<int>(H.nonZeros());'''

if old_h_block not in s:
    raise SystemExit("Could not find H construction block")

s = s.replace(old_h_block, new_h_block, 1)

# Add invalidate_derivatives and evaluate_cached_derivatives before backend_name.
anchor = "  const char* backend_name() const {"
method = '''  void invalidate_derivatives() {
    derivatives_cached_ = false;
    result_cached_ = false;
  }

  EvalResult evaluate_cached_derivatives() {
    if (!initialized_ || !derivatives_cached_ || !result_cached_) {
      throw std::runtime_error("runtime derivative cache used before initialization");
    }

    EvalResult out = cached_result_;

    if (!backend_) {
      backend_ =
          quadra::laplace::CreateLaplaceBackendForHessian(
              cached_H_,
              &recommendation_);
    }

    backend_->factorize(cached_H_);

    if (!backend_->is_spd()) {
      throw std::runtime_error("Laplace backend reported non-SPD Hessian");
    }

    out.logdet = backend_->logdet();

    const double n_x = static_cast<double>(xhat_.size());
    out.correction = 0.5 * out.logdet -
                     0.5 * n_x * std::log(2.0 * M_PI);
    out.objective = out.joint + out.correction;

    cached_result_ = out;
    result_cached_ = true;

    return out;
  }

'''

if anchor not in s:
    raise SystemExit("Could not find backend_name anchor")

s = s.replace(anchor, method + anchor, 1)

# Insert derivative-cache timing after cached_no_solve timing and before result-cache timing.
anchor_timing = '''  const auto result_cached0 = Clock::now();'''
insert_timing = '''  const auto derivative_cached0 = Clock::now();
  for (int r = 0; r < reps; ++r) {
    last = runtime.evaluate_cached_derivatives();
  }
  const auto derivative_cached1 = Clock::now();
  const double derivative_cached_total_ms =
      ms_between(derivative_cached0, derivative_cached1);
  const double derivative_cached_avg_ms =
      derivative_cached_total_ms / static_cast<double>(reps);

'''

if anchor_timing not in s:
    raise SystemExit("Could not find result-cache timing anchor")

s = s.replace(anchor_timing, insert_timing + anchor_timing, 1)

# Add derivative-cache output before result-cache output.
anchor_output = '''  std::cout << "result_cached_total_ms = " << result_cached_total_ms << "\n";'''
insert_output = '''  std::cout << "derivative_cached_total_ms = " << derivative_cached_total_ms << "\n";
  std::cout << "derivative_cached_avg_ms = " << derivative_cached_avg_ms << "\n";\n'''

if anchor_output not in s:
    raise SystemExit("Could not find result-cache output anchor")

s = s.replace(anchor_output, insert_output + anchor_output, 1)

p.write_text(s)
print("Installed derivative artifact cache benchmark.")
PY

cat <<'EOF'

Installed Phase 4 derivative artifact cache benchmark.

Run:
  ./run_state_space_surplus_latent_runtime_phase2.sh 20

Expected additional output:
  derivative_cached_total_ms = ...
  derivative_cached_avg_ms = ...

Interpretation:
  cached_avg_ms            = rebuild FD gradient/Hessian at cached xhat
  derivative_cached_avg_ms = reuse cached H/derivatives, still factorize/logdet
  result_cached_avg_ms     = return full cached EvalResult directly

EOF
