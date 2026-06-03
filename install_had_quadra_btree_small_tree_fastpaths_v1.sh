#!/usr/bin/env bash
set -euo pipefail

# install_had_quadra_btree_small_tree_fastpaths_v1.sh
#
# Low-risk BTree micro-optimization for the current generic reverse path.
#
# Motivation:
#   profile still shows hot:
#     BTree::Query
#     BTree::Insert
#     BTree::Skew/Split
#
# Hypothesis:
#   Many per-vertex sparse trees are tiny in RW1-like structures.
#
# Patch:
#   1. Add an initial reserve before first Insert so tiny trees avoid allocator
#      churn on first few nodes.
#   2. Add Query fast paths:
#        - empty tree => 0
#        - single-node tree => direct key check
#
# This keeps all tree semantics unchanged.

mkdir -p .quadra_patch_backups

target="core/had_quadra.hpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/had_quadra.hpp.btree_small_tree_fastpaths.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/btree_small_tree_fastpaths.py <<'PYEOF'
from pathlib import Path
import re

p = Path("core/had_quadra.hpp")
s = p.read_text()

# Add constants near BTree if absent.
if "QUADRA_BTREE_INITIAL_RESERVE" not in s:
    anchor = s.find("struct BTree")
    if anchor < 0:
        raise SystemExit("Could not find struct BTree")
    s = s[:anchor] + "static constexpr std::size_t QUADRA_BTREE_INITIAL_RESERVE = 4;\n\n" + s[anchor:]

# Patch BTree::Insert: add reserve on first insert.
m = re.search(r'(inline void Insert\s*\([^)]*\)\s*\{)', s)
if not m:
    raise SystemExit("Could not find BTree::Insert")

insert_start = m.end()
window = s[insert_start:insert_start + 500]
if "QUADRA_BTREE_INITIAL_RESERVE" not in window:
    reserve_code = """
    if (nodes.empty() && nodes.capacity() < QUADRA_BTREE_INITIAL_RESERVE) {
      nodes.reserve(QUADRA_BTREE_INITIAL_RESERVE);
    }
"""
    s = s[:insert_start] + reserve_code + s[insert_start:]

# Patch BTree::Query with empty/single-node fast path.
m = re.search(r'(inline Real Query\s*\([^)]*\)\s*const\s*\{)', s)
if not m:
    # Some code may use non-const Query.
    m = re.search(r'(inline Real Query\s*\([^)]*\)\s*\{)', s)

if not m:
    raise SystemExit("Could not find BTree::Query")

query_start = m.end()
window = s[query_start:query_start + 800]
if "small tree fast path" not in window:
    # Determine key parameter name from signature.
    sig = m.group(1)
    params = sig[sig.find("(")+1:sig.rfind(")")]
    # likely "const VertexId key" or "unsigned int key"
    key_name = params.strip().split()[-1].replace("&", "").replace("*", "")
    if key_name.endswith(","):
        key_name = key_name[:-1]
    if not key_name or key_name == ")":
        key_name = "key"

    fast_code = f"""
    // small tree fast path
    if (nodes.empty()) {{
      return Real(0.0);
    }}
    if (nodes.size() == 1) {{
      return nodes[0].key == {key_name} ? nodes[0].val : Real(0.0);
    }}
"""
    s = s[:query_start] + fast_code + s[query_start:]

p.write_text(s)
PYEOF

python3 /tmp/btree_small_tree_fastpaths.py

cat <<'EOF'

Installed BTree small-tree fast paths.

Run:
  ./run_had_quadra_nonzero_batch_directional_test.sh
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

Watch:
  - grad diff = 0
  - obj diff = 0
  - reverse time
  - BTree Query/Insert samples if profiling again

EOF
