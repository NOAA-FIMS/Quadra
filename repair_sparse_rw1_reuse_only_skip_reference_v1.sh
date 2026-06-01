#!/usr/bin/env bash
set -euo pipefail

# repair_sparse_rw1_reuse_only_skip_reference_v1.sh
#
# In --reuse-only mode, skip exact_gradient_rebuild() entirely.
# The earlier reuse-only mode still computed one reference result before timing,
# which polluted xctrace samples with exact_gradient_rebuild().
#
# Normal mode remains unchanged and still checks grad/objective diffs.

mkdir -p .quadra_patch_backups

target="benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp.reuse_only_skip_reference.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/reuse_only_skip_reference.py <<'PYEOF'
from pathlib import Path

p = Path("benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp")
s = p.read_text()

old = """    } else {
        // One untimed reference evaluation for correctness without polluting
        // the profile with rebuild-path samples.
        last_rebuild = exact_gradient_rebuild(
            m, K, theta, uhat, Hinv, joint_grad, joint_obj, logdet, factor);
    }"""

new = """    } else {
        // In reuse-only profiling mode, skip the rebuild reference entirely
        // so xctrace samples are not polluted by exact_gradient_rebuild().
        last_rebuild.objective = 0.0;
        last_rebuild.gradient = Eigen::VectorXd::Zero(K);
    }"""

if old not in s:
    raise SystemExit("Could not find reuse-only rebuild-reference block")

s = s.replace(old, new, 1)

old = """    out.grad_diff = (last_rebuild.gradient - last_reuse.gradient).cwiseAbs().maxCoeff();
    out.obj_diff = std::abs(last_rebuild.objective - last_reuse.objective);"""

new = """    if (reuse_only) {
        out.grad_diff = 0.0;
        out.obj_diff = 0.0;
    } else {
        out.grad_diff = (last_rebuild.gradient - last_reuse.gradient).cwiseAbs().maxCoeff();
        out.obj_diff = std::abs(last_rebuild.objective - last_reuse.objective);
    }"""

if old not in s:
    raise SystemExit("Could not find grad/objective diff block")

s = s.replace(old, new, 1)
p.write_text(s)
PYEOF

python3 /tmp/reuse_only_skip_reference.py

cat <<'EOF'

Patched --reuse-only mode to skip exact_gradient_rebuild() entirely.

Run:
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10
  ./build/benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse 10 --reuse-only

Then profile:
  xcrun xctrace record \
    --template "Time Profiler" \
    --output quadra_reuse_only.trace \
    --launch ./build/benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse \
    -- 10 --reuse-only

EOF
