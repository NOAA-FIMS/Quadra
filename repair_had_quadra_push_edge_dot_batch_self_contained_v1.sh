#!/usr/bin/env bash
set -euo pipefail

# repair_had_quadra_push_edge_dot_batch_self_contained_v1.sh
#
# Fix compile error:
#   PushEdgeDotBatch was inserted before EnsureBatchDotTreeSlot and therefore
#   cannot call it.
#
# This patch replaces that call with local bounds checks and direct access,
# making PushEdgeDotBatch self-contained.

mkdir -p .quadra_patch_backups

target="core/had_quadra.hpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/had_quadra.hpp.push_edge_dot_batch_self_contained.$(date +%Y%m%d_%H%M%S).bak"

python3 - <<'PYEOF'
from pathlib import Path

p = Path("core/had_quadra.hpp")
s = p.read_text()

old = """    BTree &tree = EnsureBatchDotTreeSlot(direction, outer, "PushEdgeDotBatch");
    tree.Insert(inner, valDot);"""

new = """    auto &trees = g_ADGraph->soEdgesDotBatch[direction];

    if (outer >= trees.size()) {
      std::cerr << "PushEdgeDotBatch tree index out of range: direction="
                << direction
                << " outer=" << outer
                << " trees.size=" << trees.size()
                << " vertices.size=" << g_ADGraph->vertices.size()
                << "\\n";
      std::abort();
    }

    trees[outer].Insert(inner, valDot);"""

if old not in s:
    raise SystemExit("Could not find EnsureBatchDotTreeSlot call in PushEdgeDotBatch")

s = s.replace(old, new, 1)
p.write_text(s)
PYEOF

cat <<'EOF'

Repaired PushEdgeDotBatch to be self-contained.

Run:
  ./run_had_quadra_nonzero_batch_directional_test.sh

EOF
