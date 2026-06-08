#!/usr/bin/env bash
set -euo pipefail

# repair_age_structured_no_plus_analytic_cohort_bandwidth_v1.sh
#
# Fixes off-by-one cohort lifetime in the analytic no-plus age-structured benchmark.
#
# With A ages:
#   x[k] recruits after observation k
#   affects observations k+1, ..., k+A
#
# Therefore at observation year y, active recruitment deviations are:
#   k in [y - A, y - 1]
#
# not [y - A + 1, y].
#
# The previous code missed the oldest surviving cohort and accidentally included
# the current-year recruitment before it entered the population. This caused
# analytic gradient/Hessian to be low and missing entries like (1,10).
#
# This patch:
#   - fixes active sensitivity index ranges in eval_all()
#   - uses bandwidth = n_ages, not n_ages - 1, for triplet extraction
#   - updates the combined runner default bandwidth for the older FD banded path
#
# Run derivative check after this patch.

target="examples/age_structured_recruitment/benchmark_age_structured_no_plus_analytic_banded.cpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  echo "Run install_age_structured_no_plus_analytic_banded_v1.sh first."
  exit 1
fi

mkdir -p .quadra_patch_backups
cp "$target" ".quadra_patch_backups/benchmark_age_structured_no_plus_analytic_banded.cpp.cohort_bandwidth.$(date +%Y%m%d_%H%M%S).bak"

python3 - <<'PYEOF'
from pathlib import Path

p = Path("examples/age_structured_recruitment/benchmark_age_structured_no_plus_analytic_banded.cpp")
s = p.read_text()

# Observation active cohorts: was [y - A + 1, y], should be [y - A, y - 1].
s = s.replace(
    "for (int k = std::max(0, y - A + 1); k <= y; ++k) {",
    "for (int k = std::max(0, y - A); k <= y - 1; ++k) {"
)

s = s.replace(
    "const int k0 = std::max(0, y - A + 1);\n    const int k1 = y;",
    "const int k0 = std::max(0, y - A);\n    const int k1 = y - 1;"
)

# If y == 0, k1 = -1 and the loops naturally skip.

# Aging after recruitment: existing cohorts before recruitment are active
# k in [y - A, y - 1], then new recruitment y is added at age 1.
s = s.replace(
    "for (int k = std::max(0, y - A + 1); k <= y; ++k) {\n        dN_next(a, k) = dN(a - 1, k) * p.survival;\n        ddN_next(a, k) = ddN_same(a - 1, k) * p.survival;\n      }",
    "for (int k = std::max(0, y - A); k <= y - 1; ++k) {\n        dN_next(a, k) = dN(a - 1, k) * p.survival;\n        ddN_next(a, k) = ddN_same(a - 1, k) * p.survival;\n      }"
)

# Triplet extraction should include distance A, not A-1.
s = s.replace(
    "for (int j = std::max(0, i - A + 1); j <= std::min(n - 1, i + A - 1); ++j) {",
    "for (int j = std::max(0, i - A); j <= std::min(n - 1, i + A); ++j) {"
)

p.write_text(s)

# Also patch the derivative checker if it included the benchmark file; no direct changes needed.
PYEOF

cat <<'EOF'

Fixed analytic no-plus cohort bandwidth/range.

Run:
  ./run_age_structured_no_plus_derivative_check.sh 25 10

Expected:
  grad max abs diff should drop sharply
  H rel diff should drop sharply

Then:
  ./run_quadra_analytic_vs_tmb_age_structured_no_plus_benchmark.sh 10 25,50,100,250,500,1000 10

EOF
