#!/usr/bin/env bash
set -euo pipefail

# install_sparse_rw1_exact_gradient_reuse_counter_breakdown_v1.sh
#
# Adds real-path batch counter reporting to:
#   benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp
#
# Counters are reset immediately before:
#   workspace.propagate_directional_batch()
#
# and accumulated immediately after it.
#
# This measures the slow exact-gradient reuse path, not the lightweight
# standalone batch counter diagnostic.

mkdir -p .quadra_patch_backups

target="benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp.real_counter_breakdown.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/quadra_reuse_counter_breakdown_patch.py <<'PYEOF'
from pathlib import Path

p = Path("benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp")
s = p.read_text()

old_row = """    double reuse_assemble_ms = 0.0;
    int vertices = 0;"""
new_row = """    double reuse_assemble_ms = 0.0;
    double avg_queries = 0.0;
    double avg_pushdots = 0.0;
    double avg_inserts = 0.0;
    int vertices = 0;"""
s = s.replace(old_row, new_row, 1)

old_sums = """    double reuse_assemble_sum = 0.0;

    const auto t2 = Clock::now();"""
new_sums = """    double reuse_assemble_sum = 0.0;
    double query_sum = 0.0;
    double pushdot_sum = 0.0;
    double insert_sum = 0.0;

    const auto t2 = Clock::now();"""
s = s.replace(old_sums, new_sums, 1)

old_reverse = """        workspace.propagate_directional_batch();
        const auto rb3 = Clock::now();"""
new_reverse = """        had::ResetBatchDirectionalCounters();
        workspace.propagate_directional_batch();
        const auto rb3 = Clock::now();

        query_sum += static_cast<double>(had::g_batch_query_count);
        pushdot_sum += static_cast<double>(had::g_batch_pushdot_count);
        insert_sum += static_cast<double>(had::g_batch_insert_count);"""
if old_reverse not in s:
    raise SystemExit("Could not find reverse timing block")
s = s.replace(old_reverse, new_reverse, 1)

old_assign = """    out.reuse_reverse_ms = reuse_reverse_sum / static_cast<double>(reps);
    out.reuse_assemble_ms = reuse_assemble_sum / static_cast<double>(reps);
    out.grad_diff = (last_rebuild.gradient - last_reuse.gradient).cwiseAbs().maxCoeff();"""
new_assign = """    out.reuse_reverse_ms = reuse_reverse_sum / static_cast<double>(reps);
    out.reuse_assemble_ms = reuse_assemble_sum / static_cast<double>(reps);
    out.avg_queries = query_sum / static_cast<double>(reps);
    out.avg_pushdots = pushdot_sum / static_cast<double>(reps);
    out.avg_inserts = insert_sum / static_cast<double>(reps);
    out.grad_diff = (last_rebuild.gradient - last_reuse.gradient).cwiseAbs().maxCoeff();"""
s = s.replace(old_assign, new_assign, 1)

old_header = """              << std::setw(14) << "assemble"
              << std::setw(16) << "grad diff"
              << std::setw(16) << "obj diff"
              << "\\n";"""
new_header = """              << std::setw(14) << "assemble"
              << std::setw(14) << "queries"
              << std::setw(14) << "pushdots"
              << std::setw(14) << "inserts"
              << std::setw(16) << "grad diff"
              << std::setw(16) << "obj diff"
              << "\\n";"""
s = s.replace(old_header, new_header, 1)

old_print = """                  << std::setw(14) << r.reuse_reverse_ms
                  << std::setw(14) << r.reuse_assemble_ms
                  << std::setw(16) << r.grad_diff
                  << std::setw(16) << r.obj_diff
                  << "\\n";"""
new_print = """                  << std::setw(14) << r.reuse_reverse_ms
                  << std::setw(14) << r.reuse_assemble_ms
                  << std::setw(14) << r.avg_queries
                  << std::setw(14) << r.avg_pushdots
                  << std::setw(14) << r.avg_inserts
                  << std::setw(16) << r.grad_diff
                  << std::setw(16) << r.obj_diff
                  << "\\n";"""
s = s.replace(old_print, new_print, 1)

p.write_text(s)
PYEOF

python3 /tmp/quadra_reuse_counter_breakdown_patch.py

cat <<'EOF'

Installed real-path reuse counter breakdown.

Run:
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

EOF
