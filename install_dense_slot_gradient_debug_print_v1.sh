#!/usr/bin/env bash
set -euo pipefail

# install_dense_slot_gradient_debug_print_v1.sh
#
# Adds a one-time debug print for m=100 comparing:
#   BTree reuse gradient
#   dense-slot gradient
#   difference
#
# Purpose:
#   dense slots are extremely fast but dense diff is large, so we need to see
#   whether the direct formula is zero, sign-flipped, missing a term, or scaled.

mkdir -p .quadra_patch_backups

target="benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp.dense_debug.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/dense_debug_patch.py <<'PYEOF'
from pathlib import Path

p = Path("benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp")
s = p.read_text()

old = """    out.dense_grad_diff =
        (last_reuse.gradient - last_dense.result.gradient).cwiseAbs().maxCoeff();

    if (reuse_only) {"""

new = """    out.dense_grad_diff =
        (last_reuse.gradient - last_dense.result.gradient).cwiseAbs().maxCoeff();

    static bool printed_dense_debug = false;
    if (!printed_dense_debug && m == 100) {
        printed_dense_debug = true;
        std::cerr << "\\nDense-slot debug m=100\\n";
        std::cerr << "reuse gradient = "
                  << last_reuse.gradient.transpose() << "\\n";
        std::cerr << "dense gradient = "
                  << last_dense.result.gradient.transpose() << "\\n";
        std::cerr << "diff gradient  = "
                  << (last_reuse.gradient - last_dense.result.gradient).transpose()
                  << "\\n";
        std::cerr << "dense traces   = "
                  << last_dense.trace_terms.transpose() << "\\n\\n";
    }

    if (reuse_only) {"""

if old not in s:
    raise SystemExit("Could not find dense diff assignment block")

s = s.replace(old, new, 1)
p.write_text(s)
PYEOF

python3 /tmp/dense_debug_patch.py

cat <<'EOF'

Installed dense-slot gradient debug print.

Run:
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 1

Paste the Dense-slot debug block.

EOF
