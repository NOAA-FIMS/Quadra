#!/usr/bin/env bash
set -euo pipefail

# repair_dense_directional_slots_match_ad_hdot_v1.sh
#
# Dense-slot debug showed exact mismatch only for:
#   log_sigma
#   log_lambda_rw
#
# reuse traces:
#   log_sigma     ~= -0.746
#   log_lambda_rw ~= -0.000435
#
# dense traces before:
#   log_sigma     ~= -103.603
#   log_lambda_rw ~= 31.3722
#
# The differences correspond to direct Huu derivative terms:
#   d/d log_sigma     exp(-2 log_sigma)
#   d/d log_lambda_rw lambda_rw RW precision terms
#
# To match the current AD Hdot path, remove those direct dense Hdot terms.
# Keep the total-u direction contribution through beta * exp(u).

mkdir -p .quadra_patch_backups

target="benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp.dense_match_ad_hdot.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/repair_dense_match_ad_hdot.py <<'PYEOF'
from pathlib import Path

p = Path("benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp")
s = p.read_text()

old_sigma = """            if (theta_dir[1] != 0.0) {
                hdot_diag += -2.0 * inv_sigma2 * theta_dir[1];
            }

"""

if old_sigma in s:
    s = s.replace(old_sigma, "", 1)
else:
    print("WARNING: direct log_sigma Hdot block not found")

old_rw_diag = """            const int rw_degree = (i > 0 ? 1 : 0) + (i + 1 < m ? 1 : 0);
            if (theta_dir[3] != 0.0) {
                hdot_diag += static_cast<double>(rw_degree) *
                             lambda_rw * theta_dir[3];
            }

"""

if old_rw_diag in s:
    s = s.replace(old_rw_diag, "", 1)
else:
    print("WARNING: direct log_lambda_rw diagonal Hdot block not found")

old_rw_sub = """                double hdot_sub = 0.0;
                if (theta_dir[3] != 0.0) {
                    hdot_sub += -lambda_rw * theta_dir[3];
                }

                trace += 2.0 * selected_inverse(i, i - 1) * hdot_sub;"""

new_rw_sub = """                const double hdot_sub = 0.0;
                trace += 2.0 * selected_inverse(i, i - 1) * hdot_sub;"""

if old_rw_sub in s:
    s = s.replace(old_rw_sub, new_rw_sub, 1)
else:
    print("WARNING: direct log_lambda_rw subdiag Hdot block not found")

p.write_text(s)
PYEOF

python3 /tmp/repair_dense_match_ad_hdot.py

cat <<'EOF'

Repaired dense directional slots to match current AD Hdot path.

Run:
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 1
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

Expected:
  dense diff near 0.
  dense ms still tiny.

EOF
