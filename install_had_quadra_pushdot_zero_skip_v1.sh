#!/usr/bin/env bash
set -euo pipefail

# install_had_quadra_pushdot_zero_skip_v1.sh
#
# Reduces PushEdgeDotBatch call volume.
#
# Current condition:
#   if (edge_dw != 0 || soDot != 0) PushEdgeDotBatch(...)
#
# But the actual value is:
#   valDot = edge_dw * soEdge.w + edge.w * soDot
#
# It can still be zero even when one input is nonzero. This patch computes
# valDot at the call site and skips PushEdgeDotBatch when valDot == 0.
#
# Correctness:
#   PushEdgeDotBatch adds valDot into Hdot storage. Skipping zero additions
#   is exact.

mkdir -p .quadra_patch_backups

target="core/had_quadra.hpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/had_quadra.hpp.pushdot_zero_skip.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/quadra_pushdot_zero_skip.py <<'PYEOF'
from pathlib import Path

p = Path("core/had_quadra.hpp")
s = p.read_text()

old1 = """        const Real e1dw_k =
            kk < e1.dwBatch.size() ? e1.dwBatch[kk] : Real(0.0);
        if (e1dw_k != Real(0.0) || soDot != Real(0.0)) {
          ++g_batch_pushdot_count;
          PushEdgeDotBatch(
              kk,
              e1,
              soEdge,
              e1dw_k,
              soDot);
        }"""

new1 = """        const Real e1dw_k =
            kk < e1.dwBatch.size() ? e1.dwBatch[kk] : Real(0.0);
        const Real e1ValDot = e1dw_k * soEdge.w + e1.w * soDot;
        if (e1ValDot != Real(0.0)) {
          ++g_batch_pushdot_count;
          PushEdgeDotBatchValue(
              kk,
              e1.to,
              soEdge.to,
              e1ValDot);
        }"""

old2 = """          const Real e2dw_k =
              kk < e2.dwBatch.size() ? e2.dwBatch[kk] : Real(0.0);
          if (e2dw_k != Real(0.0) || soDot != Real(0.0)) {
            ++g_batch_pushdot_count;
            PushEdgeDotBatch(
                kk,
                e2,
                soEdge,
                e2dw_k,
                soDot);
          }"""

new2 = """          const Real e2dw_k =
              kk < e2.dwBatch.size() ? e2.dwBatch[kk] : Real(0.0);
          const Real e2ValDot = e2dw_k * soEdge.w + e2.w * soDot;
          if (e2ValDot != Real(0.0)) {
            ++g_batch_pushdot_count;
            PushEdgeDotBatchValue(
                kk,
                e2.to,
                soEdge.to,
                e2ValDot);
          }"""

if old1 not in s:
    raise SystemExit("first PushEdgeDotBatch call block not found")
if old2 not in s:
    raise SystemExit("second PushEdgeDotBatch call block not found")

# Add PushEdgeDotBatchValue helper if not present. It avoids recomputing valDot.
if "PushEdgeDotBatchValue" not in s:
    anchor = s.find("inline void PushEdgeDotBatch(")
    if anchor < 0:
        raise SystemExit("PushEdgeDotBatch not found")

    helper = """
inline void PushEdgeDotBatchValue(const size_t direction,
                                  const VertexId fo_to,
                                  const VertexId so_to,
                                  const Real valDot) {
  if (direction >= g_ADGraph->selfSoEdgesDotBatch.size() ||
      direction >= g_ADGraph->soEdgesDotBatch.size()) {
    std::cerr << "PushEdgeDotBatchValue direction out of range: direction="
              << direction
              << " self size=" << g_ADGraph->selfSoEdgesDotBatch.size()
              << " edge size=" << g_ADGraph->soEdgesDotBatch.size()
              << "\\\\n";
    std::abort();
  }

  if (fo_to == so_to) {
    auto &selfDots = g_ADGraph->selfSoEdgesDotBatch[direction];
    if (fo_to >= selfDots.size()) {
      std::cerr << "PushEdgeDotBatchValue self index out of range: direction="
                << direction
                << " vertex=" << fo_to
                << " size=" << selfDots.size() << "\\\\n";
      std::abort();
    }

    if (!AddBatchDirectionalSlotValue(direction, fo_to, fo_to,
                                      Real(2.0) * valDot)) {
      selfDots[fo_to] += Real(2.0) * valDot;
    }
  } else {
    const VertexId outer = std::max(fo_to, so_to);
    const VertexId inner = std::min(fo_to, so_to);

    auto &trees = g_ADGraph->soEdgesDotBatch[direction];

    if (outer >= trees.size()) {
      std::cerr << "PushEdgeDotBatchValue tree index out of range: direction="
                << direction
                << " outer=" << outer
                << " trees.size=" << trees.size()
                << " vertices.size=" << g_ADGraph->vertices.size()
                << "\\\\n";
      std::abort();
    }

    if (!AddBatchDirectionalSlotValue(direction, outer, inner, valDot)) {
      trees[outer].Insert(inner, valDot);
    }
  }
}

"""
    s = s[:anchor] + helper + s[anchor:]

s = s.replace(old1, new1, 1)
s = s.replace(old2, new2, 1)

p.write_text(s)
PYEOF

python3 /tmp/quadra_pushdot_zero_skip.py

cat <<'EOF'

Installed PushEdgeDot zero-skip fast path.

Run:
  ./run_had_quadra_nonzero_batch_directional_test.sh
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

Watch:
  - grad diff remains 0
  - pushdots count should drop
  - reverse time should improve if PushEdgeDotBatch call overhead is material

EOF
