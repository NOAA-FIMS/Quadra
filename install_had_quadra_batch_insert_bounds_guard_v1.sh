#!/usr/bin/env bash
set -euo pipefail

mkdir -p .quadra_patch_backups

target="core/had_quadra.hpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/had_quadra.hpp.batch_insert_bounds_guard.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/install_batch_insert_bounds_guard.py <<'PYEOF'
from pathlib import Path

p = Path("core/had_quadra.hpp")
s = p.read_text()

if "#include <iostream>" not in s:
    if "#include <cmath>" in s:
        s = s.replace("#include <cmath>", "#include <cmath>\n#include <iostream>", 1)
    else:
        s = "#include <iostream>\n" + s

if "EnsureBatchDotTreeSlot" not in s:
    anchor = s.find("inline void PropagateAdjointDirectionalBatch()")
    if anchor < 0:
        raise SystemExit("PropagateAdjointDirectionalBatch not found")

    helper = """
inline BTree &EnsureBatchDotTreeSlot(const size_t direction,
                                     const VertexId vertex,
                                     const char *site) {
  if (direction >= g_ADGraph->soEdgesDotBatch.size()) {
    std::cerr << "soEdgesDotBatch direction out of range at "
              << site << ": direction=" << direction
              << " size=" << g_ADGraph->soEdgesDotBatch.size() << "\\n";
    std::abort();
  }

  auto &trees = g_ADGraph->soEdgesDotBatch[direction];

  if (vertex >= trees.size()) {
    std::cerr << "soEdgesDotBatch vertex out of range at "
              << site << ": direction=" << direction
              << " vertex=" << vertex
              << " trees.size=" << trees.size()
              << " vertices.size=" << g_ADGraph->vertices.size() << "\\n";
    std::abort();
  }

  return trees[vertex];
}

"""
    s = s[:anchor] + helper + s[anchor:]

replacements = [
    (
"""++g_batch_insert_count;
            g_ADGraph->soEdgesDotBatch[kk][std::max(e1.to, e2.to)].Insert(
                std::min(e1.to, e2.to), crossDot);""",
"""++g_batch_insert_count;
            EnsureBatchDotTreeSlot(
                kk, std::max(e1.to, e2.to), "crossDot")
                .Insert(std::min(e1.to, e2.to), crossDot);"""
    ),
    (
"""++g_batch_insert_count;
          g_ADGraph->soEdgesDotBatch[kk][std::max(e1.to, e2.to)].Insert(
              std::min(e1.to, e2.to), createDot);""",
"""++g_batch_insert_count;
          EnsureBatchDotTreeSlot(
              kk, std::max(e1.to, e2.to), "createDot")
              .Insert(std::min(e1.to, e2.to), createDot);"""
    ),
    (
"""g_ADGraph->soEdgesDotBatch[kk][std::max(e1.to, e2.to)].Insert(
                std::min(e1.to, e2.to), crossDot);""",
"""EnsureBatchDotTreeSlot(
                kk, std::max(e1.to, e2.to), "crossDot")
                .Insert(std::min(e1.to, e2.to), crossDot);"""
    ),
    (
"""g_ADGraph->soEdgesDotBatch[kk][std::max(e1.to, e2.to)].Insert(
              std::min(e1.to, e2.to), createDot);""",
"""EnsureBatchDotTreeSlot(
              kk, std::max(e1.to, e2.to), "createDot")
              .Insert(std::min(e1.to, e2.to), createDot);"""
    ),
]

for old, new in replacements:
    s = s.replace(old, new)

p.write_text(s)
PYEOF

python3 /tmp/install_batch_insert_bounds_guard.py

cat <<'EOF'

Installed batch insert bounds guard.

Run:
  ./run_had_quadra_nonzero_batch_directional_test.sh

If it still crashes, run:
  lldb ./build/tests/test_had_quadra_nonzero_batch_directional
  run
  bt

EOF
