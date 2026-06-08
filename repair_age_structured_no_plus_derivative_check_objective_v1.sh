#!/usr/bin/env bash
set -euo pipefail

# repair_age_structured_no_plus_derivative_check_objective_v1.sh
#
# Fix derivative checker to finite-difference eval_all(...).objective instead
# of calling a separate joint_x() function that is not present in the analytic
# benchmark source.

target="examples/age_structured_recruitment/check_age_structured_no_plus_derivatives.cpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  echo "Run install_age_structured_no_plus_derivative_check_v1.sh first."
  exit 1
fi

mkdir -p .quadra_patch_backups
cp "$target" ".quadra_patch_backups/check_age_structured_no_plus_derivatives.cpp.objective.$(date +%Y%m%d_%H%M%S).bak"

python3 - <<'PYEOF'
from pathlib import Path

p = Path("examples/age_structured_recruitment/check_age_structured_no_plus_derivatives.cpp")
s = p.read_text()

# Insert helper objective function inside anonymous namespace.
needle = "namespace {\n\n"
helper = """namespace {

double objective_x(const Data& data, const Parameters& par, const Eigen::VectorXd& x) {
  return eval_all(data, par, x).objective;
}

"""
if needle not in s:
    raise SystemExit("Could not find anonymous namespace opening")

s = s.replace(needle, helper, 1)

s = s.replace("joint_x(data, par, xp)", "objective_x(data, par, xp)")
s = s.replace("joint_x(data, par, xm)", "objective_x(data, par, xm)")
s = s.replace("joint_x(data, par, xpp)", "objective_x(data, par, xpp)")
s = s.replace("joint_x(data, par, xpm)", "objective_x(data, par, xpm)")
s = s.replace("joint_x(data, par, xmp)", "objective_x(data, par, xmp)")
s = s.replace("joint_x(data, par, xmm)", "objective_x(data, par, xmm)")
s = s.replace("const double f0 = joint_x(data, par, x);", "const double f0 = objective_x(data, par, x);")

p.write_text(s)
PYEOF

cat <<'EOF'

Repaired derivative checker to finite-difference eval_all().objective.

Run:
  ./run_age_structured_no_plus_derivative_check.sh 25 10

EOF
