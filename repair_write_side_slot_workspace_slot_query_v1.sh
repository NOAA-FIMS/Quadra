#!/usr/bin/env bash
set -euo pipefail

# repair_write_side_slot_workspace_slot_query_v1.sh
#
# Extends the experimental batch directional slot workspace so batch reverse
# reads soDot from direct slots when a mapped slot exists:
#
#   old:
#     soDot = soEdgesDotBatch[k][vid].Query(it->key)
#
#   new:
#     if mapped slot exists for (vid,it->key):
#       soDot = batchDirectionalSlotValues[k][slot]
#     else:
#       soDot = BTree.Query(...)
#
# This targets query traffic in the real reverse path.

mkdir -p .quadra_patch_backups

target="core/had_quadra.hpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/had_quadra.hpp.slot_query.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/slot_query_patch.py <<'PYEOF'
from pathlib import Path

p = Path("core/had_quadra.hpp")
s = p.read_text()

# Add helpers after GetBatchDirectionalSlotValue.
if "TryGetBatchDirectionalSlotValue" not in s:
    anchor = s.find("inline Real GetBatchDirectionalSlotValue")
    if anchor < 0:
        raise SystemExit("GetBatchDirectionalSlotValue not found")

    # insert after function end
    brace = s.find("{", anchor)
    depth = 0
    end = None
    for i in range(brace, len(s)):
        if s[i] == "{":
            depth += 1
        elif s[i] == "}":
            depth -= 1
            if depth == 0:
                end = i + 1
                break
    if end is None:
        raise SystemExit("Could not find end of GetBatchDirectionalSlotValue")

    helper = """
inline bool TryGetBatchDirectionalSlot(const VertexId i,
                                       const VertexId j,
                                       int &slot_out) {
  slot_out = -1;

  if (!g_ADGraph->useBatchDirectionalSlotWorkspace) {
    return false;
  }

  if (i == j) {
    if (i < g_ADGraph->batchSelfSlot.size()) {
      slot_out = g_ADGraph->batchSelfSlot[i];
      return slot_out >= 0;
    }
    return false;
  }

  const VertexId outer = std::max(i, j);
  const VertexId inner = std::min(i, j);

  if (outer >= g_ADGraph->batchSlotOuterInnerToSlot.size()) {
    return false;
  }

  const Real stored =
      g_ADGraph->batchSlotOuterInnerToSlot[outer].Query(inner);

  if (stored == Real(0.0)) {
    return false;
  }

  slot_out = static_cast<int>(stored) - 1;
  return slot_out >= 0;
}

inline bool TryGetBatchDirectionalSlotValue(const size_t direction,
                                            const VertexId i,
                                            const VertexId j,
                                            Real &value_out) {
  int slot = -1;
  if (!TryGetBatchDirectionalSlot(i, j, slot)) {
    return false;
  }

  if (direction >= g_ADGraph->batchDirectionalSlotValues.size()) {
    return false;
  }

  const auto &values = g_ADGraph->batchDirectionalSlotValues[direction];

  if (slot < 0 || slot >= static_cast<int>(values.size())) {
    return false;
  }

  value_out = values[static_cast<size_t>(slot)];
  return true;
}

"""
    s = s[:end] + "\n" + helper + s[end:]

old = """        BTree &btreeDot = g_ADGraph->soEdgesDotBatch[static_cast<size_t>(k)][vid];
        ++g_batch_query_count;
        const Real soDot = btreeDot.Query(it->key);"""

new = """        const size_t kk = static_cast<size_t>(k);
        Real soDot = Real(0.0);

        if (!TryGetBatchDirectionalSlotValue(kk, vid, it->key, soDot)) {
          BTree &btreeDot = g_ADGraph->soEdgesDotBatch[kk][vid];
          ++g_batch_query_count;
          soDot = btreeDot.Query(it->key);
        }"""

if old not in s:
    raise SystemExit("Could not find btreeDot Query block")

s = s.replace(old, new, 1)

# The previous active skip patch may create another 'const size_t kk' just below.
# Remove duplicate 'const size_t kk = static_cast<size_t>(k);' immediately in
# the first PushEdgeDotBatch block if present.
dup = """        const size_t kk = static_cast<size_t>(k);
        const Real e1dw_k =
            kk < e1.dwBatch.size() ? e1.dwBatch[kk] : Real(0.0);"""
dedup = """        const Real e1dw_k =
            kk < e1.dwBatch.size() ? e1.dwBatch[kk] : Real(0.0);"""
s = s.replace(dup, dedup, 1)

p.write_text(s)
PYEOF

python3 /tmp/slot_query_patch.py

cat <<'EOF'

Installed slot-backed soDot query fast path.

Run:
  ./run_had_quadra_nonzero_batch_directional_test.sh
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

Expected:
  grad diff remains 0.
  query count should drop for mapped RW1 slots.

EOF
