#!/usr/bin/env bash
set -euo pipefail

# install_dense_slot_trace_debug_print_v1.sh
#
# Extends dense-slot debug output to print:
#   joint_grad
#   reuse trace terms = 2 * (reuse gradient - joint_grad)
#   dense trace terms
#   trace differences
#
# This identifies exactly which analytic Hdot terms are missing/wrong.

mkdir -p .quadra_patch_backups

target="benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp.dense_trace_debug.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/dense_trace_debug_patch.py <<'PYEOF'
from pathlib import Path

p = Path("benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp")
s = p.read_text()

old = """        std::cerr << "dense traces   = "
                  << last_dense.trace_terms.transpose() << "\\n\\n";"""

new = """        const Eigen::VectorXd reuse_traces =
            2.0 * (last_reuse.gradient - joint_grad);

        std::cerr << "joint grad     = "
                  << joint_grad.transpose() << "\\n";
        std::cerr << "reuse traces   = "
                  << reuse_traces.transpose() << "\\n";
        std::cerr << "dense traces   = "
                  << last_dense.trace_terms.transpose() << "\\n";
        std::cerr << "trace diff     = "
                  << (reuse_traces - last_dense.trace_terms).transpose()
                  << "\\n\\n";"""

if old not in s:
    raise SystemExit("Could not find dense traces debug line")

s = s.replace(old, new, 1)
p.write_text(s)
PYEOF

python3 /tmp/dense_trace_debug_patch.py

cat <<'EOF'

Installed dense-slot trace debug print.

Run:
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 1

Paste the Dense-slot debug block again.

EOF
