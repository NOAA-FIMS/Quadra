#!/usr/bin/env bash
set -euo pipefail

# repair_flat_accumulator_backend_counter_scope_v2.sh
#
# Fixes compile errors where flat-accumulator slot hit/miss counters are used
# before they are declared.
#
# This inserts declarations immediately before AddBatchDirectionalSlotValue(),
# guaranteeing visibility for:
#   AddBatchDirectionalSlotValue()
#   TryGetBatchDirectionalSlotValue()
#
# It also avoids duplicate declarations if they already exist before that point.

mkdir -p .quadra_patch_backups

target="core/had_quadra.hpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/had_quadra.hpp.counter_scope_v2.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/repair_counter_scope_v2.py <<'PYEOF'
from pathlib import Path
import re

p = Path("core/had_quadra.hpp")
s = p.read_text()

counter_block = """// Slot-workspace hit/miss counters.
// These must appear before AddBatchDirectionalSlotValue() and
// TryGetBatchDirectionalSlotValue(), which may use them.
inline std::uint64_t g_batch_slot_query_hit_count = 0;
inline std::uint64_t g_batch_slot_query_miss_count = 0;
inline std::uint64_t g_batch_slot_write_hit_count = 0;
inline std::uint64_t g_batch_slot_write_miss_count = 0;

"""

add_pos = s.find("inline bool AddBatchDirectionalSlotValue")
if add_pos < 0:
    raise SystemExit("AddBatchDirectionalSlotValue not found")

prefix = s[:add_pos]

# If any existing declarations are after AddBatchDirectionalSlotValue, remove them
# to avoid later duplicate declarations.
s = re.sub(
    r'\n?inline std::uint64_t g_batch_slot_query_hit_count = 0;\n'
    r'inline std::uint64_t g_batch_slot_query_miss_count = 0;\n'
    r'inline std::uint64_t g_batch_slot_write_hit_count = 0;\n'
    r'inline std::uint64_t g_batch_slot_write_miss_count = 0;\n',
    '\n',
    s,
)

add_pos = s.find("inline bool AddBatchDirectionalSlotValue")
if add_pos < 0:
    raise SystemExit("AddBatchDirectionalSlotValue not found after cleanup")

# Ensure cstdint is included for uint64_t.
if "#include <cstdint>" not in s:
    if "#include <cmath>" in s:
        s = s.replace("#include <cmath>", "#include <cmath>\n#include <cstdint>", 1)
    else:
        s = "#include <cstdint>\n" + s

# Insert immediately before AddBatchDirectionalSlotValue.
s = s[:add_pos] + counter_block + s[add_pos:]

# Add reset statements only if ResetBatchDirectionalCounters exists.
reset_pos = s.find("inline void ResetBatchDirectionalCounters")
if reset_pos >= 0:
    brace = s.find("{", reset_pos)
    end = None
    depth = 0
    for i in range(brace, len(s)):
        if s[i] == "{":
            depth += 1
        elif s[i] == "}":
            depth -= 1
            if depth == 0:
                end = i
                break

    if end is not None:
        body = s[brace:end]
        if "g_batch_slot_write_miss_count = 0;" not in body:
            insert = """
  g_batch_slot_query_hit_count = 0;
  g_batch_slot_query_miss_count = 0;
  g_batch_slot_write_hit_count = 0;
  g_batch_slot_write_miss_count = 0;
"""
            s = s[:end] + insert + s[end:]

p.write_text(s)
PYEOF

python3 /tmp/repair_counter_scope_v2.py

cat <<'EOF'

Fixed flat-accumulator counter declaration scope.

Run:
  ./run_had_quadra_nonzero_batch_directional_test.sh
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

EOF
