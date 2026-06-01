#!/usr/bin/env bash
set -euo pipefail

mkdir -p .quadra_patch_backups

target="core/had_quadra.hpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/had_quadra.hpp.btree_reuse_storage.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/btree_reuse_storage_patch.py <<'PYEOF'
from pathlib import Path
import re

p = Path("core/had_quadra.hpp")
s = p.read_text()

# Patch only the first BTree Clear() block that contains nodes/root.
matches = list(re.finditer(r'inline void Clear\(\)\s*\{.*?\n\s*\}', s, re.S))
target_match = None
for m in matches:
    block = m.group(0)
    if "nodes" in block and "root" in block:
        target_match = m
        break

if target_match is None:
    raise SystemExit("Could not find BTree::Clear() block")

new_clear = """inline void Clear() {
    nodes.clear();
    root = 0;
  }"""

s = s[:target_match.start()] + new_clear + s[target_match.end():]

if "inline void Reserve(const size_t n)" not in s:
    insert = new_clear + """

  inline void Reserve(const size_t n) {
    if (nodes.capacity() < n) {
      nodes.reserve(n);
    }
  }

  inline size_t Capacity() const {
    return nodes.capacity();
  }

  inline size_t Size() const {
    return nodes.size();
  }"""
    s = s.replace(new_clear, insert, 1)

if "ReserveDirectionalBTreeStorage" not in s:
    anchor = s.find("inline void ResizeDirectionalBatch")
    if anchor < 0:
        raise SystemExit("ResizeDirectionalBatch not found")

    helper = """
inline void ReserveDirectionalBTreeStorage(const size_t reserve_per_tree) {
  for (auto &tree : g_ADGraph->soEdges) {
    tree.Reserve(reserve_per_tree);
  }
  for (auto &tree : g_ADGraph->soEdgesDot) {
    tree.Reserve(reserve_per_tree);
  }
  for (auto &trees_for_direction : g_ADGraph->soEdgesDotBatch) {
    for (auto &tree : trees_for_direction) {
      tree.Reserve(reserve_per_tree);
    }
  }
}

"""
    s = s[:anchor] + helper + s[anchor:]

p.write_text(s)
PYEOF

python3 /tmp/btree_reuse_storage_patch.py

cat <<'EOF'

Installed BTree storage reuse patch.

Run:
  ./run_had_quadra_nonzero_batch_directional_test.sh
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

EOF
