#!/usr/bin/env bash
set -euo pipefail

# repair_flat_intermediate_backend_patch_current_kk_query_v2.sh
#
# Repairs the actual current fallback query block:
#
#   BTree &btreeDot = g_ADGraph->soEdgesDotBatch[kk][vid];
#   ++g_batch_query_count;
#   soDot = btreeDot.Query(it->key);
#
# into:
#
#   if (!TryGetFlatIntermediateDirectionalValue(kk, vid, it->key, soDot)) {
#      ... BTree fallback ...
#   }
#
# This is the shape present after the previous partial flat-intermediate patch.

mkdir -p .quadra_patch_backups

target="core/had_quadra.hpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/had_quadra.hpp.current_kk_query_patch.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/patch_current_kk_query.py <<'PYEOF'
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

# Avoid double patching.
if "TryGetFlatIntermediateDirectionalValue(\n                kk, vid, static_cast<VertexId>(it->key), soDot)" in body:
    print("Current kk query block already patched.")
    p.write_text(s)
    raise SystemExit(0)

old = """          BTree &btreeDot = g_ADGraph->soEdgesDotBatch[kk][vid];
          ++g_batch_query_count;
          soDot = btreeDot.Query(it->key);"""

new = """          if (!TryGetFlatIntermediateDirectionalValue(
                  kk, vid, static_cast<VertexId>(it->key), soDot)) {
            BTree &btreeDot = g_ADGraph->soEdgesDotBatch[kk][vid];
            ++g_batch_query_count;
            soDot = btreeDot.Query(it->key);
          }"""

if old not in body:
    # More flexible fallback.
    pattern = re.compile(
        r'(?P<indent>[ \t]*)BTree\s*&\s*btreeDot\s*=\s*g_ADGraph->soEdgesDotBatch\[kk\]\[vid\]\s*;\s*\n'
        r'(?P=indent)\+\+g_batch_query_count\s*;\s*\n'
        r'(?P=indent)soDot\s*=\s*btreeDot\.Query\(it->key\)\s*;',
        re.M,
    )

    def repl(m):
        indent = m.group("indent")
        return (
            f'{indent}if (!TryGetFlatIntermediateDirectionalValue(\\n'
            f'{indent}        kk, vid, static_cast<VertexId>(it->key), soDot)) {{\\n'
            f'{indent}  BTree &btreeDot = g_ADGraph->soEdgesDotBatch[kk][vid];\\n'
            f'{indent}  ++g_batch_query_count;\\n'
            f'{indent}  soDot = btreeDot.Query(it->key);\\n'
            f'{indent}}}'
        )

    body2, n = pattern.subn(repl, body)
    if n == 0:
        candidates = []
        for line_no, line in enumerate(body.splitlines(), 1):
            if "btreeDot" in line or "soDot" in line or "soEdgesDotBatch" in line:
                candidates.append(f"{line_no}: {line}")
        raise SystemExit(
            "Could not find current kk fallback query block. Candidates:\n" +
            "\n".join(candidates[:80])
        )
    body = body2
    print(f"Patched {n} kk query block(s) by regex.")
else:
    body = body.replace(old, new, 1)
    print("Patched current kk query block by exact match.")

s = s[:start] + body + s[end:]
p.write_text(s)
PYEOF

python3 /tmp/patch_current_kk_query.py

cat <<'EOF'

Patched current kk fallback query path.

Run:
  ./run_had_quadra_nonzero_batch_directional_test.sh
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

Expected:
  - grad diff = 0
  - obj diff = 0
  - queries should drop if flat intermediate reads are used
  - reverse time should improve only if this path is dominant

EOF
