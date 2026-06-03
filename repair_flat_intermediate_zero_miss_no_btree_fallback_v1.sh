#!/usr/bin/env bash
set -euo pipefail

# repair_flat_intermediate_zero_miss_no_btree_fallback_v1.sh
#
# Optimization after diagnostics:
#
#   flat intermediate: read_hit=27500 read_miss=14995 write_hit=59990 write_miss=0
#
# Interpretation:
#   - all writes now go to the flat intermediate backend
#   - a read miss means no flat slot/value exists for that edge
#   - because write_miss=0, there should be no corresponding BTree dot value
#   - therefore the missed read can safely remain zero without BTree fallback
#
# This patch changes the current query fallback block from:
#
#   if (!TryGetFlatIntermediateDirectionalValue(..., soDot)) {
#       BTree query fallback
#   }
#
# to:
#
#   if (flat backend enabled) {
#       TryGetFlatIntermediateDirectionalValue(..., soDot);
#       // miss leaves soDot = 0
#   } else {
#       BTree query fallback
#   }
#
# Expected:
#   - queries drop from ~14995 to ~0
#   - grad diff / obj diff remain 0
#   - reverse time improves further

mkdir -p .quadra_patch_backups

target="core/had_quadra.hpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/had_quadra.hpp.zero_miss_no_btree.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/zero_miss_no_btree.py <<'PYEOF'
from pathlib import Path
import re

p = Path("core/had_quadra.hpp")
s = p.read_text()

def find_function(src, name):
    idx = src.find(name)
    if idx < 0:
        raise SystemExit(f"function {name} not found")
    start = src.rfind("\n", 0, idx) + 1
    brace = src.find("{", idx)
    depth = 0
    end = None
    for i in range(brace, len(src)):
        if src[i] == "{":
            depth += 1
        elif src[i] == "}":
            depth -= 1
            if depth == 0:
                end = i + 1
                break
    if end is None:
        raise SystemExit(f"could not find end of {name}")
    return start, end

start, end = find_function(s, "PropagateAdjointDirectionalBatch")
body = s[start:end]

old = """          if (!TryGetFlatIntermediateDirectionalValue(
                  kk, vid, static_cast<VertexId>(it->key), soDot)) {
            BTree &btreeDot = g_ADGraph->soEdgesDotBatch[kk][vid];
            ++g_batch_query_count;
            soDot = btreeDot.Query(it->key);
          }"""

new = """          if (g_ADGraph->useFlatIntermediateDirectionalBackend) {
            // A flat miss means no directional value was written for this edge.
            // Since all patched writes use the flat backend, leave soDot = 0.
            (void)TryGetFlatIntermediateDirectionalValue(
                kk, vid, static_cast<VertexId>(it->key), soDot);
          } else {
            BTree &btreeDot = g_ADGraph->soEdgesDotBatch[kk][vid];
            ++g_batch_query_count;
            soDot = btreeDot.Query(it->key);
          }"""

if old not in body:
    candidates = []
    for line_no, line in enumerate(body.splitlines(), 1):
        if "TryGetFlatIntermediateDirectionalValue" in line or "btreeDot" in line or "soDot = btreeDot.Query" in line:
            candidates.append(f"{line_no}: {line}")
    raise SystemExit("Could not find zero-miss fallback block. Candidates:\n" + "\n".join(candidates[:80]))

body = body.replace(old, new, 1)
s = s[:start] + body + s[end:]

p.write_text(s)
PYEOF

python3 /tmp/zero_miss_no_btree.py

cat <<'EOF'

Installed zero-miss no-BTree fallback optimization.

Run:
  ./run_had_quadra_nonzero_batch_directional_test.sh
  ./run_flat_intermediate_backend_verification_suite.sh

Or minimally:
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

Watch:
  - queries should drop
  - grad diff = 0
  - obj diff = 0
  - reverse time

EOF
