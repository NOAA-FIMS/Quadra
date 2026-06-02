#!/usr/bin/env bash
set -euo pipefail

mkdir -p .quadra_patch_backups

target="core/had_quadra.hpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/had_quadra.hpp.batch_push_edge_dot.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/repair_batch_push_edge_dot.py <<'PYEOF'
from pathlib import Path

p = Path("core/had_quadra.hpp")
s = p.read_text()

if "PushEdgeDotBatch" not in s:
    anchor = s.find("inline Real GetAdjointDotBatch")
    if anchor < 0:
        anchor = s.find("inline BTree &EnsureBatchDotTreeSlot")
    if anchor < 0:
        raise SystemExit("Could not find insertion anchor for PushEdgeDotBatch")

    helper = """
inline void PushEdgeDotBatch(const size_t direction,
                             const ADEdge &foEdge,
                             const ADEdge &soEdge,
                             const Real foEdgeDot,
                             const Real soEdgeDot) {
  const Real valDot = foEdgeDot * soEdge.w + foEdge.w * soEdgeDot;

  if (direction >= g_ADGraph->selfSoEdgesDotBatch.size() ||
      direction >= g_ADGraph->soEdgesDotBatch.size()) {
    std::cerr << "PushEdgeDotBatch direction out of range: direction="
              << direction
              << " self size=" << g_ADGraph->selfSoEdgesDotBatch.size()
              << " edge size=" << g_ADGraph->soEdgesDotBatch.size()
              << "\\n";
    std::abort();
  }

  if (foEdge.to == soEdge.to) {
    auto &selfDots = g_ADGraph->selfSoEdgesDotBatch[direction];
    if (foEdge.to >= selfDots.size()) {
      std::cerr << "PushEdgeDotBatch self index out of range: direction="
                << direction
                << " vertex=" << foEdge.to
                << " size=" << selfDots.size() << "\\n";
      std::abort();
    }
    selfDots[foEdge.to] += Real(2.0) * valDot;
  } else {
    const VertexId outer = std::max(foEdge.to, soEdge.to);
    const VertexId inner = std::min(foEdge.to, soEdge.to);

    BTree &tree = EnsureBatchDotTreeSlot(direction, outer, "PushEdgeDotBatch");
    tree.Insert(inner, valDot);
  }
}

"""
    s = s[:anchor] + helper + s[anchor:]

old1 = """        // Reuse scalar PushEdgeDot by temporarily assigning dw.
        const Real oldDw1 = e1.dw;
        e1.dw = e1.dwBatch[static_cast<size_t>(k)];
        ++g_batch_pushdot_count;
        PushEdgeDot(e1, soEdge, soDot);
        e1.dw = oldDw1;"""

new1 = """        ++g_batch_pushdot_count;
        PushEdgeDotBatch(
            static_cast<size_t>(k),
            e1,
            soEdge,
            e1.dwBatch[static_cast<size_t>(k)],
            soDot);"""

if old1 in s:
    s = s.replace(old1, new1, 1)
else:
    print("WARNING: first scalar PushEdgeDot batch block not found")

old2 = """          const Real oldDw2 = e2.dw;
          e2.dw = e2.dwBatch[static_cast<size_t>(k)];
          ++g_batch_pushdot_count;
          PushEdgeDot(e2, soEdge, soDot);
          e2.dw = oldDw2;"""

new2 = """          ++g_batch_pushdot_count;
          PushEdgeDotBatch(
              static_cast<size_t>(k),
              e2,
              soEdge,
              e2.dwBatch[static_cast<size_t>(k)],
              soDot);"""

if old2 in s:
    s = s.replace(old2, new2, 1)
else:
    print("WARNING: second scalar PushEdgeDot batch block not found")

p.write_text(s)
PYEOF

python3 /tmp/repair_batch_push_edge_dot.py

cat <<'EOF'

Installed batch-specific PushEdgeDotBatch repair.

Run:
  ./run_had_quadra_nonzero_batch_directional_test.sh
  ./run_had_quadra_directional_batch_propagation_test.sh

EOF
