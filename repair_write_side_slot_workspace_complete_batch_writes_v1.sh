#!/usr/bin/env bash
set -euo pipefail

# repair_write_side_slot_workspace_complete_batch_writes_v1.sh
#
# The first slot-workspace patch only redirected PushEdgeDotBatch() writes.
# The batch reverse also writes Hdot terms directly:
#
#   selfSoEdgesDotBatch[kk][e1.to] += ...
#   selfSoEdgesDotBatch[kk][e2.to] += ...
#   selfSoEdgesDotBatch[kk][e1.to] += 2*crossDot
#   selfSoEdgesDotBatch[kk][e1.to] += createDot
#   selfSoEdgesDotBatch[kk][e1.to] += 2*createDot
#   soEdgesDotBatch[...] Insert(crossDot/createDot)
#
# When exact_gradient() reads only slots, missing these direct writes causes
# the big grad diff. This patch routes those direct writes through
# AddBatchDirectionalSlotValue() with fallback to the old storage.

mkdir -p .quadra_patch_backups

target="core/had_quadra.hpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/had_quadra.hpp.complete_batch_slot_writes.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/complete_batch_slot_writes.py <<'PYEOF'
from pathlib import Path

p = Path("core/had_quadra.hpp")
s = p.read_text()

repls = [
(
"""        g_ADGraph->selfSoEdgesDotBatch[kk][e1.to] +=
            Real(2.0) * e1.w * e1.dwBatch[kk] * S + e1.w * e1.w * SDot;""",
"""        const Real e1SelfDot =
            Real(2.0) * e1.w * e1.dwBatch[kk] * S + e1.w * e1.w * SDot;
        if (!AddBatchDirectionalSlotValue(kk, e1.to, e1.to, e1SelfDot)) {
          g_ADGraph->selfSoEdgesDotBatch[kk][e1.to] += e1SelfDot;
        }"""
),
(
"""          g_ADGraph->selfSoEdgesDotBatch[kk][e2.to] +=
              Real(2.0) * e2.w * e2.dwBatch[kk] * S + e2.w * e2.w * SDot;""",
"""          const Real e2SelfDot =
              Real(2.0) * e2.w * e2.dwBatch[kk] * S + e2.w * e2.w * SDot;
          if (!AddBatchDirectionalSlotValue(kk, e2.to, e2.to, e2SelfDot)) {
            g_ADGraph->selfSoEdgesDotBatch[kk][e2.to] += e2SelfDot;
          }"""
),
(
"""            g_ADGraph->selfSoEdgesDotBatch[kk][e1.to] += Real(2.0) * crossDot;""",
"""            if (!AddBatchDirectionalSlotValue(
                    kk, e1.to, e1.to, Real(2.0) * crossDot)) {
              g_ADGraph->selfSoEdgesDotBatch[kk][e1.to] += Real(2.0) * crossDot;
            }"""
),
(
"""          g_ADGraph->selfSoEdgesDotBatch[kk][e1.to] += createDot;""",
"""          if (!AddBatchDirectionalSlotValue(kk, e1.to, e1.to, createDot)) {
            g_ADGraph->selfSoEdgesDotBatch[kk][e1.to] += createDot;
          }"""
),
(
"""          g_ADGraph->selfSoEdgesDotBatch[kk][e1.to] += Real(2.0) * createDot;""",
"""          if (!AddBatchDirectionalSlotValue(
                  kk, e1.to, e1.to, Real(2.0) * createDot)) {
            g_ADGraph->selfSoEdgesDotBatch[kk][e1.to] += Real(2.0) * createDot;
          }"""
),
]

changed = 0
for old, new in repls:
    if old in s:
        s = s.replace(old, new)
        changed += 1
    else:
        print("WARNING: pattern not found:", old.splitlines()[0])

# Ensure explicit offdiag createDot block has guarded slot writes.
# Existing crossDot/createDot offdiag writes may already call EnsureBatchDotTreeSlot.
# Replace those with AddBatchDirectionalSlotValue-first guarded path.
old = """            ++g_batch_insert_count;
            EnsureBatchDotTreeSlot(
                kk, std::max(e1.to, e2.to), "crossDot")
                .Insert(std::min(e1.to, e2.to), crossDot);"""
new = """            if (!AddBatchDirectionalSlotValue(kk, e1.to, e2.to, crossDot)) {
              ++g_batch_insert_count;
              EnsureBatchDotTreeSlot(
                  kk, std::max(e1.to, e2.to), "crossDot")
                  .Insert(std::min(e1.to, e2.to), crossDot);
            }"""
if old in s:
    s = s.replace(old, new)
    changed += 1
else:
    print("WARNING: crossDot offdiag pattern not found")

old = """          ++g_batch_insert_count;
          EnsureBatchDotTreeSlot(
              kk, std::max(e1.to, e2.to), "createDot")
              .Insert(std::min(e1.to, e2.to), createDot);"""
new = """          if (!AddBatchDirectionalSlotValue(kk, e1.to, e2.to, createDot)) {
            ++g_batch_insert_count;
            EnsureBatchDotTreeSlot(
                kk, std::max(e1.to, e2.to), "createDot")
                .Insert(std::min(e1.to, e2.to), createDot);
          }"""
if old in s:
    s = s.replace(old, new)
    changed += 1
else:
    print("WARNING: createDot offdiag pattern not found")

if changed == 0:
    raise SystemExit("No batch write patterns were changed")

p.write_text(s)
print(f"changed patterns: {changed}")
PYEOF

python3 /tmp/complete_batch_slot_writes.py

cat <<'EOF'

Completed write-side slot routing for direct batch Hdot writes.

Run:
  ./run_had_quadra_nonzero_batch_directional_test.sh
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

Expected:
  grad diff should return to 0.
  insert count should drop if mapped writes are captured by slots.

EOF
