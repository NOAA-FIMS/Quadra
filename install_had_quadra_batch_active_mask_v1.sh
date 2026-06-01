#!/usr/bin/env bash
set -euo pipefail

# install_had_quadra_batch_active_mask_v1.sh
#
# Adds a cached per-vertex active direction mask for batched directional reverse.
#
# Current hotspot:
#   BatchDirectionHasLocalSignal(vertex,e1,e2,vid,k)
#
# is called many times inside PropagateAdjointDirectionalBatch().
#
# This patch:
#   - adds ADVertex::batchActiveDirectionMask
#   - computes the mask after forward replay / storage normalization
#   - replaces repeated BatchDirectionHasLocalSignal(...) checks with:
#
#       if ((vertex.batchActiveDirectionMask & (1ULL << kk)) == 0) continue;
#
# Conservative:
#   - supports up to 64 directions through the fast mask
#   - falls back to BatchDirectionHasLocalSignal for directions >= 64
#   - keeps the existing helper available

mkdir -p .quadra_patch_backups

target="core/had_quadra.hpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/had_quadra.hpp.batch_active_mask.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/quadra_batch_active_mask_patch.py <<'PYEOF'
from pathlib import Path

p = Path("core/had_quadra.hpp")
s = p.read_text()

# Add include for uint64_t if not present.
if "#include <cstdint>" not in s:
    if "#include <cmath>" in s:
        s = s.replace("#include <cmath>", "#include <cmath>\n#include <cstdint>", 1)
    else:
        s = "#include <cstdint>\n" + s

# Add field to ADVertex.
if "batchActiveDirectionMask" not in s:
    old = """        // Batched directional derivative of soW along seeded primal tangents.
        std::vector<Real> soWDotBatch;

  // Replay primal value and operation metadata."""
    new = """        // Batched directional derivative of soW along seeded primal tangents.
        std::vector<Real> soWDotBatch;

        // Cached active-direction mask for batched directional reverse.
        // Bit k is set when direction k has local signal at this vertex.
        // Used for fast paths when nBatchDirections <= 64.
        std::uint64_t batchActiveDirectionMask = 0;

  // Replay primal value and operation metadata."""
    if old not in s:
        raise SystemExit("Could not find ADVertex soWDotBatch block")
    s = s.replace(old, new, 1)

# Reset in ResizeDirectionalBatch/ClearDirectionalBatch loops.
old = """            v.e2.dwBatch.assign(static_cast<size_t>(nDirections), Real(0.0));"""
new = """            v.e2.dwBatch.assign(static_cast<size_t>(nDirections), Real(0.0));
            v.batchActiveDirectionMask = 0;"""
s = s.replace(old, new)

# Add helper functions after BatchDirectionHasLocalSignal.
if "ComputeBatchActiveDirectionMasks" not in s:
    anchor = s.find("inline bool BatchDirectionHasLocalSignal")
    if anchor < 0:
        raise SystemExit("BatchDirectionHasLocalSignal not found")

    # find function end
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
        raise SystemExit("Could not find end of BatchDirectionHasLocalSignal")

    helper = """
inline void ComputeBatchActiveDirectionMasks(const int nDirections) {
  const int cappedDirections = std::min(nDirections, 64);

  for (VertexId vid = 0;
       vid < static_cast<VertexId>(g_ADGraph->vertices.size());
       ++vid) {
    ADVertex &vertex = g_ADGraph->vertices[vid];
    ADEdge &e1 = vertex.e1;
    ADEdge &e2 = vertex.e2;

    std::uint64_t mask = 0;

    for (int k = 0; k < cappedDirections; ++k) {
      const size_t kk = static_cast<size_t>(k);
      if (BatchDirectionHasLocalSignal(vertex, e1, e2, vid, kk)) {
        mask |= (std::uint64_t(1) << kk);
      }
    }

    vertex.batchActiveDirectionMask = mask;
  }
}

inline bool BatchDirectionMaskHasSignal(const ADVertex &vertex,
                                        const ADVertex & /*unused*/,
                                        const ADEdge &e1,
                                        const ADEdge &e2,
                                        const VertexId vid,
                                        const size_t kk) {
  if (kk < 64) {
    return (vertex.batchActiveDirectionMask & (std::uint64_t(1) << kk)) != 0;
  }

  return BatchDirectionHasLocalSignal(vertex, e1, e2, vid, kk);
}

"""
    s = s[:end] + "\n" + helper + s[end:]

# Call ComputeBatchActiveDirectionMasks after replay/storage normalization before reverse sweep.
marker = """  for (VertexId vid = n_vertices - 1; vid > 0; --vid) {"""
idx = s.find("inline void PropagateAdjointDirectionalBatch()")
pos = s.find(marker, idx)
if pos < 0:
    raise SystemExit("Could not find reverse sweep loop")
if "ComputeBatchActiveDirectionMasks(nDirections);" not in s[idx:pos]:
    insert = """  ComputeBatchActiveDirectionMasks(nDirections);

"""
    s = s[:pos] + insert + s[pos:]

# Replace exact checks in reverse sweep.
old = """      if (!BatchDirectionHasLocalSignal(vertex, e1, e2, vid, kk)) {
        continue;
      }"""
new = """      if (!BatchDirectionMaskHasSignal(vertex, vertex, e1, e2, vid, kk)) {
        continue;
      }"""
s = s.replace(old, new)

p.write_text(s)
PYEOF

python3 /tmp/quadra_batch_active_mask_patch.py

cat <<'EOF'

Installed HAD Quadra batch active-direction mask v1.

Run:
  ./run_had_quadra_nonzero_batch_directional_test.sh
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

Watch:
  - grad diff remains 0
  - reverse time changes

EOF
