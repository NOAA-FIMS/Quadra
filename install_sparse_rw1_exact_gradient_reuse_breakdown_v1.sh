#!/usr/bin/env bash
set -euo pipefail

# install_sparse_rw1_exact_gradient_reuse_breakdown_v1.sh
#
# Adds reuse-path timing buckets to:
#   benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp
#
# Buckets:
#   base adjoint
#   seed directions
#   batch reverse
#   assemble/trace
#
# Purpose:
#   The standalone Hdot pipeline says batch reverse is small, but the
#   exact-gradient reuse path is still ~30 ms at m=500. This benchmark
#   identifies which reuse sub-step is responsible.

mkdir -p .quadra_patch_backups

target="benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp.reuse_breakdown.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/quadra_reuse_breakdown_patch.py <<'PYEOF'
from pathlib import Path

p = Path("benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp")
s = p.read_text()

# Add fields to Row.
old_row = """struct Row {
    double rebuild_ms = 0.0;
    double reuse_ms = 0.0;
    double speedup = 0.0;
    double grad_diff = 0.0;
    double obj_diff = 0.0;
    int vertices = 0;
};"""

new_row = """struct Row {
    double rebuild_ms = 0.0;
    double reuse_ms = 0.0;
    double speedup = 0.0;
    double grad_diff = 0.0;
    double obj_diff = 0.0;
    double reuse_base_ms = 0.0;
    double reuse_seed_ms = 0.0;
    double reuse_reverse_ms = 0.0;
    double reuse_assemble_ms = 0.0;
    int vertices = 0;
};"""

if old_row in s:
    s = s.replace(old_row, new_row, 1)

# Replace reuse timing loop.
old_loop = """    const auto t2 = Clock::now();
    for (int r = 0; r < reps; ++r) {
        workspace.propagate_base_adjoint();
        workspace.seed_directions(theta, uhat, factor);
        workspace.propagate_directional_batch();

        last_reuse = workspace.exact_gradient(
            selected_inverse, joint_grad, joint_obj, logdet);
    }
    const auto t3 = Clock::now();"""

new_loop = """    double reuse_base_sum = 0.0;
    double reuse_seed_sum = 0.0;
    double reuse_reverse_sum = 0.0;
    double reuse_assemble_sum = 0.0;

    const auto t2 = Clock::now();
    for (int r = 0; r < reps; ++r) {
        const auto rb0 = Clock::now();
        workspace.propagate_base_adjoint();
        const auto rb1 = Clock::now();

        workspace.seed_directions(theta, uhat, factor);
        const auto rb2 = Clock::now();

        workspace.propagate_directional_batch();
        const auto rb3 = Clock::now();

        last_reuse = workspace.exact_gradient(
            selected_inverse, joint_grad, joint_obj, logdet);
        const auto rb4 = Clock::now();

        reuse_base_sum += ms_between(rb0, rb1);
        reuse_seed_sum += ms_between(rb1, rb2);
        reuse_reverse_sum += ms_between(rb2, rb3);
        reuse_assemble_sum += ms_between(rb3, rb4);
    }
    const auto t3 = Clock::now();"""

if old_loop not in s:
    raise SystemExit("Could not find reuse timing loop to replace")
s = s.replace(old_loop, new_loop, 1)

# Add assignment after reuse_ms/speedup.
old_assign = """    out.rebuild_ms = ms_between(t0, t1) / static_cast<double>(reps);
    out.reuse_ms = ms_between(t2, t3) / static_cast<double>(reps);
    out.speedup = out.rebuild_ms / out.reuse_ms;
    out.grad_diff = (last_rebuild.gradient - last_reuse.gradient).cwiseAbs().maxCoeff();
    out.obj_diff = std::abs(last_rebuild.objective - last_reuse.objective);"""

new_assign = """    out.rebuild_ms = ms_between(t0, t1) / static_cast<double>(reps);
    out.reuse_ms = ms_between(t2, t3) / static_cast<double>(reps);
    out.speedup = out.rebuild_ms / out.reuse_ms;
    out.reuse_base_ms = reuse_base_sum / static_cast<double>(reps);
    out.reuse_seed_ms = reuse_seed_sum / static_cast<double>(reps);
    out.reuse_reverse_ms = reuse_reverse_sum / static_cast<double>(reps);
    out.reuse_assemble_ms = reuse_assemble_sum / static_cast<double>(reps);
    out.grad_diff = (last_rebuild.gradient - last_reuse.gradient).cwiseAbs().maxCoeff();
    out.obj_diff = std::abs(last_rebuild.objective - last_reuse.objective);"""

if old_assign not in s:
    raise SystemExit("Could not find assignment block")
s = s.replace(old_assign, new_assign, 1)

# Update header print.
old_header = """    std::cout << std::setw(8) << "m"
              << std::setw(12) << "vertices"
              << std::setw(16) << "rebuild ms"
              << std::setw(16) << "reuse ms"
              << std::setw(14) << "speedup"
              << std::setw(16) << "grad diff"
              << std::setw(16) << "obj diff"
              << "\\n";"""

new_header = """    std::cout << std::setw(8) << "m"
              << std::setw(12) << "vertices"
              << std::setw(16) << "rebuild ms"
              << std::setw(16) << "reuse ms"
              << std::setw(14) << "speedup"
              << std::setw(14) << "base"
              << std::setw(14) << "seed"
              << std::setw(14) << "reverse"
              << std::setw(14) << "assemble"
              << std::setw(16) << "grad diff"
              << std::setw(16) << "obj diff"
              << "\\n";"""

if old_header in s:
    s = s.replace(old_header, new_header, 1)

old_print = """        std::cout << std::setw(8) << m
                  << std::setw(12) << r.vertices
                  << std::setw(16) << r.rebuild_ms
                  << std::setw(16) << r.reuse_ms
                  << std::setw(14) << r.speedup
                  << std::setw(16) << r.grad_diff
                  << std::setw(16) << r.obj_diff
                  << "\\n";"""

new_print = """        std::cout << std::setw(8) << m
                  << std::setw(12) << r.vertices
                  << std::setw(16) << r.rebuild_ms
                  << std::setw(16) << r.reuse_ms
                  << std::setw(14) << r.speedup
                  << std::setw(14) << r.reuse_base_ms
                  << std::setw(14) << r.reuse_seed_ms
                  << std::setw(14) << r.reuse_reverse_ms
                  << std::setw(14) << r.reuse_assemble_ms
                  << std::setw(16) << r.grad_diff
                  << std::setw(16) << r.obj_diff
                  << "\\n";"""

if old_print in s:
    s = s.replace(old_print, new_print, 1)

p.write_text(s)
PYEOF

python3 /tmp/quadra_reuse_breakdown_patch.py

cat <<'EOF'

Installed sparse RW1 exact-gradient reuse breakdown.

Patched:
  benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp

Run:
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

EOF
