#!/usr/bin/env bash
set -euo pipefail

# install_batch_slot_hit_miss_counters_v1.sh
#
# Adds counters for the experimental slot-query fast path:
#
#   g_batch_slot_query_hit_count
#   g_batch_slot_query_miss_count
#   g_batch_slot_write_hit_count
#   g_batch_slot_write_miss_count
#
# and prints them in benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp.
#
# Purpose:
#   The query count did not drop after adding slot-backed soDot reads.
#   This will tell us whether the slot map is missing because reverse sweep
#   operates on intermediate AD vertices rather than final random-effect vertices.

mkdir -p .quadra_patch_backups

header="core/had_quadra.hpp"
bench="benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp"

if [[ ! -f "$header" ]]; then
  echo "ERROR: missing $header"
  exit 1
fi
if [[ ! -f "$bench" ]]; then
  echo "ERROR: missing $bench"
  exit 1
fi

cp "$header" ".quadra_patch_backups/had_quadra.hpp.slot_hit_miss_counters.$(date +%Y%m%d_%H%M%S).bak"
cp "$bench" ".quadra_patch_backups/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp.slot_hit_miss_counters.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/slot_hit_miss_patch.py <<'PYEOF'
from pathlib import Path

h = Path("core/had_quadra.hpp")
s = h.read_text()

if "g_batch_slot_query_hit_count" not in s:
    old = """inline std::uint64_t g_batch_insert_count = 0;"""
    new = """inline std::uint64_t g_batch_insert_count = 0;
inline std::uint64_t g_batch_slot_query_hit_count = 0;
inline std::uint64_t g_batch_slot_query_miss_count = 0;
inline std::uint64_t g_batch_slot_write_hit_count = 0;
inline std::uint64_t g_batch_slot_write_miss_count = 0;"""
    s = s.replace(old, new, 1)

old = """  g_batch_insert_count = 0;"""
new = """  g_batch_insert_count = 0;
  g_batch_slot_query_hit_count = 0;
  g_batch_slot_query_miss_count = 0;
  g_batch_slot_write_hit_count = 0;
  g_batch_slot_write_miss_count = 0;"""
s = s.replace(old, new, 1)

# Count write hit/miss in AddBatchDirectionalSlotValue.
old = """  if (!g_ADGraph->useBatchDirectionalSlotWorkspace) {
    return false;
  }"""
new = """  if (!g_ADGraph->useBatchDirectionalSlotWorkspace) {
    ++g_batch_slot_write_miss_count;
    return false;
  }"""
s = s.replace(old, new, 1)

old = """  if (direction >= g_ADGraph->batchDirectionalSlotValues.size()) {
    return false;
  }"""
new = """  if (direction >= g_ADGraph->batchDirectionalSlotValues.size()) {
    ++g_batch_slot_write_miss_count;
    return false;
  }"""
s = s.replace(old, new, 1)

old = """  if (slot < 0) {
    return false;
  }"""
new = """  if (slot < 0) {
    ++g_batch_slot_write_miss_count;
    return false;
  }"""
s = s.replace(old, new, 1)

old = """  if (slot >= static_cast<int>(values.size())) {
    return false;
  }

  values[static_cast<size_t>(slot)] += value;
  return true;"""
new = """  if (slot >= static_cast<int>(values.size())) {
    ++g_batch_slot_write_miss_count;
    return false;
  }

  values[static_cast<size_t>(slot)] += value;
  ++g_batch_slot_write_hit_count;
  return true;"""
s = s.replace(old, new, 1)

# Count query hit/miss in TryGetBatchDirectionalSlotValue.
old = """  int slot = -1;
  if (!TryGetBatchDirectionalSlot(i, j, slot)) {
    return false;
  }"""
new = """  int slot = -1;
  if (!TryGetBatchDirectionalSlot(i, j, slot)) {
    ++g_batch_slot_query_miss_count;
    return false;
  }"""
s = s.replace(old, new, 1)

old = """  if (direction >= g_ADGraph->batchDirectionalSlotValues.size()) {
    return false;
  }"""
new = """  if (direction >= g_ADGraph->batchDirectionalSlotValues.size()) {
    ++g_batch_slot_query_miss_count;
    return false;
  }"""
s = s.replace(old, new, 1)

old = """  if (slot < 0 || slot >= static_cast<int>(values.size())) {
    return false;
  }

  value_out = values[static_cast<size_t>(slot)];
  return true;"""
new = """  if (slot < 0 || slot >= static_cast<int>(values.size())) {
    ++g_batch_slot_query_miss_count;
    return false;
  }

  value_out = values[static_cast<size_t>(slot)];
  ++g_batch_slot_query_hit_count;
  return true;"""
s = s.replace(old, new, 1)

h.write_text(s)

b = Path("benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp")
t = b.read_text()

# Add fields.
old = """    double avg_inserts = 0.0;
    int vertices = 0;"""
new = """    double avg_inserts = 0.0;
    double avg_slot_q_hit = 0.0;
    double avg_slot_q_miss = 0.0;
    double avg_slot_w_hit = 0.0;
    double avg_slot_w_miss = 0.0;
    int vertices = 0;"""
t = t.replace(old, new, 1)

old = """    double insert_sum = 0.0;

    const auto t2 = Clock::now();"""
new = """    double insert_sum = 0.0;
    double slot_q_hit_sum = 0.0;
    double slot_q_miss_sum = 0.0;
    double slot_w_hit_sum = 0.0;
    double slot_w_miss_sum = 0.0;

    const auto t2 = Clock::now();"""
t = t.replace(old, new, 1)

old = """        insert_sum += static_cast<double>(had::g_batch_insert_count);"""
new = """        insert_sum += static_cast<double>(had::g_batch_insert_count);
        slot_q_hit_sum += static_cast<double>(had::g_batch_slot_query_hit_count);
        slot_q_miss_sum += static_cast<double>(had::g_batch_slot_query_miss_count);
        slot_w_hit_sum += static_cast<double>(had::g_batch_slot_write_hit_count);
        slot_w_miss_sum += static_cast<double>(had::g_batch_slot_write_miss_count);"""
t = t.replace(old, new, 1)

old = """    out.avg_inserts = insert_sum / static_cast<double>(reps);
    out.grad_diff = (last_rebuild.gradient - last_reuse.gradient).cwiseAbs().maxCoeff();"""
new = """    out.avg_inserts = insert_sum / static_cast<double>(reps);
    out.avg_slot_q_hit = slot_q_hit_sum / static_cast<double>(reps);
    out.avg_slot_q_miss = slot_q_miss_sum / static_cast<double>(reps);
    out.avg_slot_w_hit = slot_w_hit_sum / static_cast<double>(reps);
    out.avg_slot_w_miss = slot_w_miss_sum / static_cast<double>(reps);
    out.grad_diff = (last_rebuild.gradient - last_reuse.gradient).cwiseAbs().maxCoeff();"""
t = t.replace(old, new, 1)

old = """              << std::setw(14) << "inserts"
              << std::setw(16) << "grad diff\""" 
new = """              << std::setw(14) << "inserts"
              << std::setw(14) << "qhit"
              << std::setw(14) << "qmiss"
              << std::setw(14) << "whit"
              << std::setw(14) << "wmiss"
              << std::setw(16) << "grad diff\""" 
t = t.replace(old, new, 1)

old = """                  << std::setw(14) << r.avg_inserts
                  << std::setw(16) << r.grad_diff"""
new = """                  << std::setw(14) << r.avg_inserts
                  << std::setw(14) << r.avg_slot_q_hit
                  << std::setw(14) << r.avg_slot_q_miss
                  << std::setw(14) << r.avg_slot_w_hit
                  << std::setw(14) << r.avg_slot_w_miss
                  << std::setw(16) << r.grad_diff"""
t = t.replace(old, new, 1)

b.write_text(t)
PYEOF

python3 /tmp/slot_hit_miss_patch.py

cat <<'EOF'

Installed slot hit/miss counters.

Run:
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

EOF
