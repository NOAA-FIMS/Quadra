#!/usr/bin/env bash
set -euo pipefail

# install_state_space_latent_runtime_result_cache_v1.sh
#
# Adds a result cache to the Phase 3 runtime example.
#
# New mode:
#   evaluate_cached_result()
#
# This returns the last EvalResult directly when theta/xhat are unchanged.
# It demonstrates the upper bound for fixed-theta repeated evaluations:
#   no optimizer
#   no FD gradient
#   no FD Hessian
#   no backend factorization
#
# This is intentionally for the fixed-theta runtime benchmark only.
# A production evaluator should invalidate this result cache when theta changes.

target="examples/state_space_surplus_production/laplace_state_space_surplus_latent_runtime.cpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

mkdir -p .quadra_patch_backups
cp "$target" ".quadra_patch_backups/laplace_state_space_surplus_latent_runtime.cpp.result_cache.$(date +%Y%m%d_%H%M%S).bak"

python3 - "$target" <<'PY'
import sys
from pathlib import Path

p = Path(sys.argv[1])
s = p.read_text()

if "evaluate_cached_result()" in s:
    print("Result cache patch already appears to be installed.")
    raise SystemExit(0)

# Add fields to EvalResult so default cached result is safe to initialize.
# No struct change needed; EvalResult already has defaults.

# Add result cache members.
old_members = '''  bool initialized_;

  quadra::laplace::BackendRecommendation recommendation_;
  std::unique_ptr<quadra::laplace::LaplaceBackend> backend_;'''

new_members = '''  bool initialized_;
  bool result_cached_ = false;
  EvalResult cached_result_;

  quadra::laplace::BackendRecommendation recommendation_;
  std::unique_ptr<quadra::laplace::LaplaceBackend> backend_;'''

if old_members not in s:
    raise SystemExit("Could not find runtime member block")

s = s.replace(old_members, new_members, 1)

# Store result in evaluate_at_xhat before return.
old_return = '''    out.objective = out.joint + out.correction;

    return out;
  }'''

new_return = '''    out.objective = out.joint + out.correction;

    cached_result_ = out;
    result_cached_ = true;

    return out;
  }'''

if old_return not in s:
    raise SystemExit("Could not find evaluate_at_xhat return block")

s = s.replace(old_return, new_return, 1)

# Add evaluate_cached_result after evaluate_cached_no_solve().
old_cached_no_solve = '''  EvalResult evaluate_cached_no_solve() {
    if (!initialized_) {
      throw std::runtime_error("runtime cache used before initialization");
    }
    return evaluate_at_xhat();
  }

  const char* backend_name() const {'''

new_cached_no_solve = '''  EvalResult evaluate_cached_no_solve() {
    if (!initialized_) {
      throw std::runtime_error("runtime cache used before initialization");
    }
    return evaluate_at_xhat();
  }

  EvalResult evaluate_cached_result() const {
    if (!initialized_ || !result_cached_) {
      throw std::runtime_error("runtime result cache used before initialization");
    }
    return cached_result_;
  }

  const char* backend_name() const {'''

if old_cached_no_solve not in s:
    raise SystemExit("Could not find evaluate_cached_no_solve block")

s = s.replace(old_cached_no_solve, new_cached_no_solve, 1)

# Insert benchmark timing block after cached_no_solve block in main.
old_main_cached = '''  const auto cached0 = Clock::now();
  for (int r = 0; r < reps; ++r) {
    last = runtime.evaluate_cached_no_solve();
  }
  const auto cached1 = Clock::now();
  const double cached_total_ms = ms_between(cached0, cached1);
  const double cached_avg_ms = cached_total_ms / static_cast<double>(reps);

  std::cout << std::fixed << std::setprecision(6);'''

new_main_cached = '''  const auto cached0 = Clock::now();
  for (int r = 0; r < reps; ++r) {
    last = runtime.evaluate_cached_no_solve();
  }
  const auto cached1 = Clock::now();
  const double cached_total_ms = ms_between(cached0, cached1);
  const double cached_avg_ms = cached_total_ms / static_cast<double>(reps);

  const auto result_cached0 = Clock::now();
  for (int r = 0; r < reps; ++r) {
    last = runtime.evaluate_cached_result();
  }
  const auto result_cached1 = Clock::now();
  const double result_cached_total_ms =
      ms_between(result_cached0, result_cached1);
  const double result_cached_avg_ms =
      result_cached_total_ms / static_cast<double>(reps);

  std::cout << std::fixed << std::setprecision(6);'''

if old_main_cached not in s:
    raise SystemExit("Could not find main cached timing block")

s = s.replace(old_main_cached, new_main_cached, 1)

# Add output lines.
old_output = '''  std::cout << "cached_total_ms = " << cached_total_ms << "\n";
  std::cout << "cached_avg_ms = " << cached_avg_ms << "\n";'''

new_output = '''  std::cout << "cached_total_ms = " << cached_total_ms << "\n";
  std::cout << "cached_avg_ms = " << cached_avg_ms << "\n";
  std::cout << "result_cached_total_ms = " << result_cached_total_ms << "\n";
  std::cout << "result_cached_avg_ms = " << result_cached_avg_ms << "\n";'''

if old_output not in s:
    raise SystemExit("Could not find cached output block")

s = s.replace(old_output, new_output, 1)

p.write_text(s)
PY

cat <<'EOF'

Installed runtime result cache benchmark.

Run:
  ./run_state_space_surplus_latent_runtime_phase2.sh 20

Expected:
  backend = tridiagonal
  objective = -10.642176

Interpretation:
  warm_avg_ms           = runtime warm path for unchanged theta
  cached_avg_ms         = recompute FD gradient/Hessian/backend at cached xhat
  result_cached_avg_ms  = return stored result directly, valid only when theta/xhat unchanged

EOF
