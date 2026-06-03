#!/usr/bin/env bash
set -euo pipefail

# repair_flat_intermediate_backend_patch_actual_query_v1.sh
#
# The first flat intermediate backend installed, but reported:
#   WARNING: did not find primary query pattern
#
# Result:
#   correctness held, but the main soEdgesDotBatch Query path stayed on BTree,
#   so reverse time did not improve.
#
# This repair:
#   1. Removes temporary intermediate registry print from the benchmark.
#   2. Patches the actual current soEdgesDotBatch[...] Query pattern by regex.
#   3. Leaves fallback BTree behavior if a flat slot is missing.

mkdir -p .quadra_patch_backups

header="core/had_quadra.hpp"
bench="benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp"

if [[ ! -f "$header" ]]; then
  echo "ERROR: missing $header"
  exit 1
fi
if [[ ! -f "$bench" ]]; then
  echo "ERROR: missing $bench"
  exit 1
fi

cp "$header" ".quadra_patch_backups/had_quadra.hpp.flat_query_repair.$(date +%Y%m%d_%H%M%S).bak"
cp "$bench" ".quadra_patch_backups/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp.remove_registry_print.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/repair_flat_query.py <<'PYEOF'
from pathlib import Path
import re

# Remove temporary benchmark print.
b = Path("benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp")
t = b.read_text()
t = t.replace("""        if (m == 500 && r == reps - 1) {
            had::PrintIntermediateEdgeSlotRegistryDiagnostic();
        }

""", "")
b.write_text(t)

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

# Match the current query pattern flexibly:
#
#   BTree &btreeDot = g_ADGraph->soEdgesDotBatch[...][vid];
#   ++g_batch_query_count;
#   const Real soDot = btreeDot.Query(it->key);
#
# or:
#
#   auto &btreeDot = ...
#   const Real soDot = ...
#
pattern = re.compile(
    r'(?P<indent>[ \t]*)'
    r'(?:BTree|auto)\s*&\s*(?P<tree>\w+)\s*=\s*'
    r'g_ADGraph->soEdgesDotBatch\[(?P<dir>[^\]]+)\]\[vid\]\s*;\s*\n'
    r'(?P=indent)\+\+g_batch_query_count\s*;\s*\n'
    r'(?P=indent)const\s+Real\s+(?P<var>\w+)\s*=\s*(?P=tree)\.Query\(it->key\)\s*;',
    re.M,
)

def repl(m):
    indent = m.group("indent")
    tree = m.group("tree")
    direction_expr = m.group("dir")
    var = m.group("var")
    return (
        f'{indent}Real {var} = Real(0.0);\n'
        f'{indent}if (!TryGetFlatIntermediateDirectionalValue(\n'
        f'{indent}        {direction_expr}, vid, static_cast<VertexId>(it->key), {var})) {{\n'
        f'{indent}  BTree &{tree} = g_ADGraph->soEdgesDotBatch[{direction_expr}][vid];\n'
        f'{indent}  ++g_batch_query_count;\n'
        f'{indent}  {var} = {tree}.Query(it->key);\n'
        f'{indent}}}'
    )

body2, n = pattern.subn(repl, body)

if n == 0:
    # Dump likely candidates to help debugging, but fail so we do not pretend.
    candidates = []
    for line_no, line in enumerate(body.splitlines(), 1):
        if "soEdgesDotBatch" in line or ".Query(it->key)" in line:
            candidates.append(f"{line_no}: {line}")
    msg = "Could not patch query pattern. Candidates:\\n" + "\\n".join(candidates[:40])
    raise SystemExit(msg)

s = s[:start] + body2 + s[end:]
p.write_text(s)

print(f"Patched {n} soEdgesDotBatch Query pattern(s).")
PYEOF

python3 /tmp/repair_flat_query.py

cat <<'EOF'

Repaired flat intermediate backend query path.

Run:
  ./run_had_quadra_nonzero_batch_directional_test.sh
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

Watch:
  - grad diff = 0
  - obj diff = 0
  - queries should drop if flat reads are being used
  - reverse time should improve if this path is hot

EOF
