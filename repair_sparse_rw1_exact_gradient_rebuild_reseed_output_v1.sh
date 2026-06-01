#!/usr/bin/env bash
set -euo pipefail

# repair_sparse_rw1_exact_gradient_rebuild_reseed_output_v1.sh
#
# After HadGraphWorkspace started reseeding the output adjoint before batch
# reverse, the reuse path now computes nonzero trace terms. The rebuild
# reference path in benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp still
# calls had::PropagateAdjointDirectionalBatch() directly without reseeding
# output.w, so it remains the old zero-Hdot reference.
#
# This patches the rebuild benchmark path to reseed:
#
#   graph.vertices[f.varId].w = 1.0;
#
# immediately before had::PropagateAdjointDirectionalBatch().

mkdir -p .quadra_patch_backups

target="benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp.rebuild_reseed_output.$(date +%Y%m%d_%H%M%S).bak"

python3 - <<'PYEOF'
from pathlib import Path

p = Path("benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp")
s = p.read_text()

needle = """    had::PropagateAdjointDirectionalBatch();

    return compute_exact_gradient_from_current_graph("""
replacement = """    // PropagateAdjoint() consumes/clears reverse adjoints while building the
    // base Hessian. Match HadGraphWorkspace behavior by reseeding the output
    // adjoint before the directional reverse sweep.
    graph.vertices[f.varId].w = Real(1.0);
    had::PropagateAdjointDirectionalBatch();

    return compute_exact_gradient_from_current_graph("""

if needle not in s:
    raise SystemExit("Could not find rebuild direct PropagateAdjointDirectionalBatch call")

s = s.replace(needle, replacement, 1)
p.write_text(s)
PYEOF

cat <<'EOF'

Patched sparse RW1 exact-gradient rebuild reference to reseed output adjoint.

Run:
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

Expected:
  grad diff returns near 0, obj diff remains 0.
  Rebuild time may increase because it is now doing real Hdot work.

EOF
