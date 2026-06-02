#!/usr/bin/env bash
set -euo pipefail

# repair_had_quadra_batch_dot_storage_resize_v2.sh

mkdir -p .quadra_patch_backups

target="core/had_quadra.hpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/had_quadra.hpp.batch_dot_storage_resize_v2.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/repair_batch_dot_storage_resize_v2.py <<'PYEOF'
from pathlib import Path

p = Path("core/had_quadra.hpp")
s = p.read_text()

idx = s.find("inline void PropagateAdjointDirectionalBatch()")
if idx < 0:
    raise SystemExit("PropagateAdjointDirectionalBatch not found")

marker = """  PropagateDirectionalBatchForwardReplay();

  for (VertexId vid = n_vertices - 1; vid > 0; --vid) {"""

replacement = """  PropagateDirectionalBatchForwardReplay();

  // Defensive storage normalization before reverse sweep. The directional
  // reverse pass inserts into soEdgesDotBatch[direction][vertex], so every
  // direction must own a full vector of BTree slots.
  g_ADGraph->soEdgesDotBatch.resize(static_cast<size_t>(nDirections));
  g_ADGraph->selfSoEdgesDotBatch.resize(static_cast<size_t>(nDirections));

  for (int k = 0; k < nDirections; ++k) {
    auto &edgeDots = g_ADGraph->soEdgesDotBatch[static_cast<size_t>(k)];
    if (edgeDots.size() < g_ADGraph->vertices.size()) {
      edgeDots.resize(g_ADGraph->vertices.size());
    }

    auto &selfDots = g_ADGraph->selfSoEdgesDotBatch[static_cast<size_t>(k)];
    if (selfDots.size() < g_ADGraph->vertices.size()) {
      selfDots.resize(g_ADGraph->vertices.size(), Real(0.0));
    }
  }

  for (VertexId vid = n_vertices - 1; vid > 0; --vid) {"""

pos = s.find(marker, idx)
if pos < 0:
    # Fallback: insert immediately before the reverse loop marker.
    loop = "  for (VertexId vid = n_vertices - 1; vid > 0; --vid) {"
    pos = s.find(loop, idx)
    if pos < 0:
        raise SystemExit("Could not find reverse vid loop marker")

    insert = """  // Defensive storage normalization before reverse sweep. The directional
  // reverse pass inserts into soEdgesDotBatch[direction][vertex], so every
  // direction must own a full vector of BTree slots.
  g_ADGraph->soEdgesDotBatch.resize(static_cast<size_t>(nDirections));
  g_ADGraph->selfSoEdgesDotBatch.resize(static_cast<size_t>(nDirections));

  for (int k = 0; k < nDirections; ++k) {
    auto &edgeDots = g_ADGraph->soEdgesDotBatch[static_cast<size_t>(k)];
    if (edgeDots.size() < g_ADGraph->vertices.size()) {
      edgeDots.resize(g_ADGraph->vertices.size());
    }

    auto &selfDots = g_ADGraph->selfSoEdgesDotBatch[static_cast<size_t>(k)];
    if (selfDots.size() < g_ADGraph->vertices.size()) {
      selfDots.resize(g_ADGraph->vertices.size(), Real(0.0));
    }
  }

"""
    s = s[:pos] + insert + s[pos:]
else:
    s = s.replace(marker, replacement, 1)

p.write_text(s)
PYEOF

python3 /tmp/repair_batch_dot_storage_resize_v2.py

cat <<'EOF'

Repaired batch dot storage sizing before reverse sweep.

Run:
  ./run_had_quadra_nonzero_batch_directional_test.sh

EOF
