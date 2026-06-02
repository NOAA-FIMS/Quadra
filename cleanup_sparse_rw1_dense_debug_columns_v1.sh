#!/usr/bin/env bash
set -euo pipefail

# cleanup_sparse_rw1_dense_debug_columns_v1.sh
#
# Cleans benchmark-only dense validation/debug output from:
#   benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp
#
# Keeps:
#   core/laplace/structure_aware_rw1_hdot.hpp
#   compute_rw1_dense_slot_exact_gradient()
#
# Rationale:
#   The structure-aware helper is now validated:
#     dense diff = 0
#     grad diff  = 0
#     obj diff   = 0
#
# This cleanup removes noisy benchmark columns/debug printing while preserving
# the reusable helper and correctness path internally.

mkdir -p .quadra_patch_backups

target="benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp.cleanup_dense_debug.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/cleanup_dense_debug_columns.py <<'PYEOF'
from pathlib import Path
import re

p = Path("benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp")
s = p.read_text()

# Remove Row dense columns.
s = s.replace("    double dense_slots_ms = 0.0;\n", "")
s = s.replace("    double dense_speedup_vs_reuse = 0.0;\n", "")
s = s.replace("    double dense_grad_diff = 0.0;\n", "")

# Keep last_dense calculation for internal validation? No, remove timed dense loop
# and validation assignments from benchmark output. The reusable helper remains
# in the file and can be used by dedicated tests later.
s = s.replace("    DenseSlotExactGradientResult last_dense;\n", "")

# Remove dense timing loop.
s = re.sub(
    r'\n    const auto td0 = Clock::now\(\);\n'
    r'    for \(int r = 0; r < reps; \+\+r\) \{\n'
    r'        last_dense = compute_rw1_dense_slot_exact_gradient\(\n'
    r'            m, K, theta, uhat, selected_inverse, joint_grad,\n'
    r'            joint_obj, logdet, factor\);\n'
    r'    \}\n'
    r'    const auto td1 = Clock::now\(\);\n',
    '\n',
    s,
)

# Remove dense assignment block.
s = re.sub(
    r'\n    out\.dense_slots_ms = ms_between\(td0, td1\) / static_cast<double>\(reps\);\n'
    r'    out\.dense_speedup_vs_reuse =\n'
    r'        out\.dense_slots_ms > 0\.0 \? out\.reuse_ms / out\.dense_slots_ms : 0\.0;\n'
    r'    out\.dense_grad_diff =\n'
    r'        \(last_reuse\.gradient - last_dense\.result\.gradient\)\.cwiseAbs\(\)\.maxCoeff\(\);\n',
    '\n',
    s,
)

# Remove dense debug block.
s = re.sub(
    r'\n    static bool printed_dense_debug = false;\n'
    r'    if \(!printed_dense_debug && m == 100\) \{\n'
    r'.*?'
    r'    \}\n\n'
    r'    if \(reuse_only\) \{',
    '\n    if (reuse_only) {',
    s,
    flags=re.S,
)

# Restore reuse-only grad/objective diff to no-reference behavior.
s = re.sub(
    r'    if \(reuse_only\) \{\n'
    r'        out\.grad_diff = out\.dense_grad_diff;\n'
    r'        out\.obj_diff = std::abs\(last_reuse\.objective - last_dense\.result\.objective\);\n'
    r'    \} else \{',
    '    if (reuse_only) {\n'
    '        out.grad_diff = 0.0;\n'
    '        out.obj_diff = 0.0;\n'
    '    } else {',
    s,
)

# Remove dense header columns.
s = s.replace('              << std::setw(14) << "dense ms"\n', '')
s = s.replace('              << std::setw(14) << "dense spd"\n', '')
s = s.replace('              << std::setw(16) << "dense diff"\n', '')

# Remove dense print columns.
s = s.replace("                  << std::setw(14) << r.dense_slots_ms\n", "")
s = s.replace("                  << std::setw(14) << r.dense_speedup_vs_reuse\n", "")
s = s.replace("                  << std::setw(16) << r.dense_grad_diff\n", "")

p.write_text(s)
PYEOF

python3 /tmp/cleanup_dense_debug_columns.py

cat <<'EOF'

Cleaned dense validation/debug columns from sparse RW1 graph reuse benchmark.

Run:
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10
  ./build/benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse 10 --reuse-only

Expected:
  - no Dense-slot debug block
  - no dense ms/spd/diff columns
  - grad diff = 0, obj diff = 0 in normal mode
  - reuse-only still runs for profiling

EOF
