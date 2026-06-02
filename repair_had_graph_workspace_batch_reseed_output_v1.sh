#!/usr/bin/env bash
set -euo pipefail

# repair_had_graph_workspace_batch_reseed_output_v1.sh
#
# Low-level nonzero diagnostic now passes after manually reseeding output.w = 1
# before PropagateAdjointDirectionalBatch().
#
# This patch moves that behavior into HadGraphWorkspace:
#
#   PropagateAdjointDirectionalBatch()
#     Activate()
#     graph_.vertices[output_var_id_].w = 1
#     had::PropagateAdjointDirectionalBatch()
#
# This should allow ExactGradientWorkspace/LaplaceExactGradientEvaluator traces
# to become nonzero without each caller manually reseeding.

mkdir -p .quadra_patch_backups

target="core/had_graph_workspace.hpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/had_graph_workspace.hpp.batch_reseed_output.$(date +%Y%m%d_%H%M%S).bak"

python3 - <<'PYEOF'
from pathlib import Path

p = Path("core/had_graph_workspace.hpp")
s = p.read_text()

old = """    void PropagateAdjointDirectionalBatch() {
        Activate();
        had::PropagateAdjointDirectionalBatch();
    }"""

new = """    void PropagateAdjointDirectionalBatch() {
        if (!built_) {
            throw std::logic_error(
                "HadGraphWorkspace::PropagateAdjointDirectionalBatch called before Build.");
        }

        Activate();

        if (output_var_id_ >= graph_.vertices.size()) {
            throw std::out_of_range(
                "HadGraphWorkspace::PropagateAdjointDirectionalBatch output_var_id out of range.");
        }

        // PropagateAdjoint() consumes/clears reverse adjoints while building
        // the base Hessian. The directional reverse sweep needs the output
        // adjoint seed restored.
        graph_.vertices[output_var_id_].w = had::Real(1.0);

        had::PropagateAdjointDirectionalBatch();
    }"""

if old not in s:
    raise SystemExit("Could not find HadGraphWorkspace::PropagateAdjointDirectionalBatch block")

s = s.replace(old, new, 1)
p.write_text(s)
PYEOF

cat <<'EOF'

Patched HadGraphWorkspace batch reverse to reseed output adjoint.

Run:
  ./run_had_graph_workspace_test.sh
  ./run_exact_gradient_workspace_test.sh
  ./run_laplace_exact_gradient_evaluator_test.sh
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

EOF
