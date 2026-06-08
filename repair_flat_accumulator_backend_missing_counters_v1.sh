#!/usr/bin/env bash
set -euo pipefail

# repair_flat_accumulator_backend_missing_counters_v1.sh
#
# The flat-accumulator backend v2 references slot hit/miss counters that may
# not exist if the earlier diagnostics patch was not applied:
#
#   g_batch_slot_write_miss_count
#   g_batch_slot_write_hit_count
#   g_batch_slot_query_miss_count
#   g_batch_slot_query_hit_count
#
# This repair adds those counters near the existing batch directional counters
# and resets them in ResetBatchDirectionalCounters().

mkdir -p .quadra_patch_backups

target="core/had_quadra.hpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/had_quadra.hpp.missing_flat_counters.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/repair_flat_counters.py <<'PYEOF'
from pathlib import Path

p = Path("core/had_quadra.hpp")
s = p.read_text()

# Add missing counter declarations after existing batch counters.
if "g_batch_slot_write_miss_count" not in s:
    anchor = "inline std::uint64_t g_batch_insert_count = 0;"
    if anchor not in s:
        raise SystemExit("Could not find g_batch_insert_count declaration")

    s = s.replace(
        anchor,
        anchor + """
inline std::uint64_t g_batch_slot_query_hit_count = 0;
inline std::uint64_t g_batch_slot_query_miss_count = 0;
inline std::uint64_t g_batch_slot_write_hit_count = 0;
inline std::uint64_t g_batch_slot_write_miss_count = 0;""",
        1,
    )

# Reset counters inside ResetBatchDirectionalCounters().
reset_anchor = """  g_batch_insert_count = 0;"""
if reset_anchor not in s:
    raise SystemExit("Could not find ResetBatchDirectionalCounters insert reset")

reset_block = """  g_batch_insert_count = 0;
  g_batch_slot_query_hit_count = 0;
  g_batch_slot_query_miss_count = 0;
  g_batch_slot_write_hit_count = 0;
  g_batch_slot_write_miss_count = 0;"""

if "g_batch_slot_write_miss_count = 0;" not in s:
    s = s.replace(reset_anchor, reset_block, 1)

p.write_text(s)
PYEOF

python3 /tmp/repair_flat_counters.py

cat <<'EOF'

Added missing flat-accumulator slot hit/miss counters.

Run:
  ./run_had_quadra_nonzero_batch_directional_test.sh
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10
  ./run_sparse_rw1_flat_accumulator_trace_check.sh 10

EOF
