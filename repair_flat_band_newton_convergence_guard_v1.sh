#!/usr/bin/env bash
set -euo pipefail

# repair_flat_band_newton_convergence_guard_v1.sh
# Stops Newton before line search when gradient/step/objective are already converged.

files=(
  "examples/age_structured_recruitment/benchmark_age_structured_no_plus_flat_band.cpp"
  "examples/age_structured_recruitment/benchmark_age_structured_no_plus_flat_band_newton_diagnostics.cpp"
  "examples/age_structured_recruitment/benchmark_age_structured_no_plus_flat_band_alpha_trace.cpp"
  "examples/age_structured_recruitment/benchmark_age_structured_no_plus_flat_band_warmstart.cpp"
)

mkdir -p .quadra_patch_backups

for target in "${files[@]}"; do
  [[ -f "$target" ]] || continue
  cp "$target" ".quadra_patch_backups/$(basename "$target").convergence_guard.$(date +%Y%m%d_%H%M%S).bak"
  python3 - "$target" <<'PY'
import sys
from pathlib import Path

p = Path(sys.argv[1])
s = p.read_text()

# Loosen ultra-strict gradient tolerance. The benchmark prints 0.000000,
# but the old 1e-8 guard was letting near-machine-noise iterations continue.
s = s.replace('if (gnorm < 1e-8) break;', 'if (gnorm < 1e-6) break;')

# Add step-norm guard right after step_norm is available.
old = '''    const Eigen::VectorXd dx = ldlt.solve(e.gradient);
    const double step_norm = dx.norm();
    if (ldlt.info() != Eigen::Success || !dx.allFinite()) {'''
new = '''    const Eigen::VectorXd dx = ldlt.solve(e.gradient);
    const double step_norm = dx.norm();
    if (step_norm < 1e-10 * (1.0 + x.norm())) {
      break;
    }
    if (ldlt.info() != Eigen::Success || !dx.allFinite()) {'''
if old in s:
    s = s.replace(old, new, 1)

# If this file does not have step_norm instrumentation yet, add it.
old2 = '''    const Eigen::VectorXd dx = ldlt.solve(e.gradient);
    if (ldlt.info() != Eigen::Success || !dx.allFinite()) {'''
new2 = '''    const Eigen::VectorXd dx = ldlt.solve(e.gradient);
    const double step_norm = dx.norm();
    if (step_norm < 1e-10 * (1.0 + x.norm())) {
      break;
    }
    if (ldlt.info() != Eigen::Success || !dx.allFinite()) {'''
if old2 in s:
    s = s.replace(old2, new2, 1)

# Permit no strict improvement if objective is numerically unchanged and gradient is tiny.
old3 = '''      if (std::isfinite(f_candidate) && f_candidate < f) {'''
new3 = '''      const double rel_tol = 1e-12 * (1.0 + std::abs(f));
      if (std::isfinite(f_candidate) &&
          (f_candidate < f || (std::abs(f_candidate - f) <= rel_tol && gnorm < 1e-6))) {'''
if old3 in s:
    s = s.replace(old3, new3, 1)

# In alpha trace files current_f may be used; keep compatibility.
p.write_text(s)
PY
done

cat <<'EOF'

Installed Newton convergence guard.

Run diagnostics again:
  ./run_quadra_age_structured_no_plus_flat_band_newton_diagnostics.sh 3 25,50,100,250,500,1000 10

Then benchmark:
  ./run_quadra_flat_band_vs_tmb_age_structured_no_plus_benchmark.sh 10 25,50,100,250,500,1000 10

EOF
