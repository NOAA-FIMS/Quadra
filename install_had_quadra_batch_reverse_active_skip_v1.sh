#!/usr/bin/env bash
set -euo pipefail

mkdir -p .quadra_patch_backups

target="core/had_quadra.hpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/had_quadra.hpp.batch_reverse_active_skip.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/quadra_batch_reverse_active_skip.py <<'PYEOF'
from pathlib import Path

p = Path("core/had_quadra.hpp")
s = p.read_text()

if "BatchDirectionHasLocalSignal" not in s:
    anchor = s.find("inline void PropagateAdjointDirectionalBatch()")
    if anchor < 0:
        raise SystemExit("PropagateAdjointDirectionalBatch not found")

    helper = """
inline bool BatchDirectionHasLocalSignal(const ADVertex &vertex,
                                         const ADEdge &e1,
                                         const ADEdge &e2,
                                         const VertexId vid,
                                         const size_t kk) {
  const Real e1dw =
      kk < e1.dwBatch.size() ? e1.dwBatch[kk] : Real(0.0);
  const Real e2dw =
      kk < e2.dwBatch.size() ? e2.dwBatch[kk] : Real(0.0);
  const Real aDot =
      kk < vertex.wDotBatch.size() ? vertex.wDotBatch[kk] : Real(0.0);
  const Real soWDot =
      kk < vertex.soWDotBatch.size() ? vertex.soWDotBatch[kk] : Real(0.0);

  Real selfDot = Real(0.0);
  if (kk < g_ADGraph->selfSoEdgesDotBatch.size() &&
      vid < g_ADGraph->selfSoEdgesDotBatch[kk].size()) {
    selfDot = g_ADGraph->selfSoEdgesDotBatch[kk][vid];
  }

  return e1dw != Real(0.0) ||
         (e2.to != vid && e2dw != Real(0.0)) ||
         aDot != Real(0.0) ||
         soWDot != Real(0.0) ||
         selfDot != Real(0.0);
}

"""
    s = s[:anchor] + helper + s[anchor:]

old1 = """        ++g_batch_pushdot_count;
        PushEdgeDotBatch(
            static_cast<size_t>(k),
            e1,
            soEdge,
            e1.dwBatch[static_cast<size_t>(k)],
            soDot);"""
new1 = """        const size_t kk = static_cast<size_t>(k);
        const Real e1dw_k =
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
s = s.replace(old1, new1, 1)

old2 = """          ++g_batch_pushdot_count;
          PushEdgeDotBatch(
              static_cast<size_t>(k),
              e2,
              soEdge,
              e2.dwBatch[static_cast<size_t>(k)],
              soDot);"""
new2 = """          const Real e2dw_k =
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
s = s.replace(old2, new2, 1)

patterns = [
(
"""    for (int k = 0; k < nDirections; ++k) {
      const size_t kk = static_cast<size_t>(k);
      const Real SDot = g_ADGraph->selfSoEdgesDotBatch[kk][vid];

      if (S != Real(0.0) || SDot != Real(0.0)) {""",
"""    for (int k = 0; k < nDirections; ++k) {
      const size_t kk = static_cast<size_t>(k);
      if (!BatchDirectionHasLocalSignal(vertex, e1, e2, vid, kk)) {
        continue;
      }

      const Real SDot = g_ADGraph->selfSoEdgesDotBatch[kk][vid];

      if (S != Real(0.0) || SDot != Real(0.0)) {"""
),
(
"""    for (int k = 0; k < nDirections; ++k) {
      const size_t kk = static_cast<size_t>(k);
      const Real aDot = vertex.wDotBatch[kk];

      if ((a != Real(0.0) || aDot != Real(0.0)) &&
          (vertex.soW != Real(0.0) || vertex.soWDotBatch[kk] != Real(0.0))) {""",
"""    for (int k = 0; k < nDirections; ++k) {
      const size_t kk = static_cast<size_t>(k);
      if (!BatchDirectionHasLocalSignal(vertex, e1, e2, vid, kk)) {
        continue;
      }

      const Real aDot = vertex.wDotBatch[kk];

      if ((a != Real(0.0) || aDot != Real(0.0)) &&
          (vertex.soW != Real(0.0) || vertex.soWDotBatch[kk] != Real(0.0))) {"""
),
(
"""    for (int k = 0; k < nDirections; ++k) {
      const size_t kk = static_cast<size_t>(k);
      const Real aDot = vertex.wDotBatch[kk];

      if (a != Real(0.0) || aDot != Real(0.0)) {""",
"""    for (int k = 0; k < nDirections; ++k) {
      const size_t kk = static_cast<size_t>(k);
      if (!BatchDirectionHasLocalSignal(vertex, e1, e2, vid, kk)) {
        continue;
      }

      const Real aDot = vertex.wDotBatch[kk];

      if (a != Real(0.0) || aDot != Real(0.0)) {"""
),
]

for old, new in patterns:
    if old in s:
        s = s.replace(old, new, 1)
    else:
        print("WARNING: pattern not found")

p.write_text(s)
PYEOF

python3 /tmp/quadra_batch_reverse_active_skip.py

cat <<'EOF'

Installed HAD Quadra batch reverse active-direction skip v1.

Run:
  ./run_had_quadra_nonzero_batch_directional_test.sh
  ./run_had_quadra_directional_batch_propagation_test.sh
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

EOF
