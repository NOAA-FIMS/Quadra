#!/usr/bin/env bash
set -euo pipefail

mkdir -p .quadra_patch_backups

header="core/had_quadra.hpp"
accum="core/had/batch_directional_flat_accumulator.hpp"

if [[ ! -f "$header" ]]; then
  echo "ERROR: missing $header"
  exit 1
fi

if [[ ! -f "$accum" ]]; then
  echo "ERROR: missing $accum"
  exit 1
fi

cp "$header" ".quadra_patch_backups/had_quadra.hpp.preallocate_flat_slots.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/preallocate_flat_slots.py <<'PYEOF'
from pathlib import Path

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

if "EstimateFlatIntermediateDirectionalSlotCapacity" not in s:
    anchor = s.find("inline void EnableFlatIntermediateDirectionalBackend")
    if anchor < 0:
        raise SystemExit("EnableFlatIntermediateDirectionalBackend not found")

    helper = '''
inline std::size_t EstimateFlatIntermediateDirectionalSlotCapacity() {
  std::size_t edge_count = 0;

  for (const auto &tree : g_ADGraph->soEdges) {
    edge_count += tree.Size();
  }

  const std::size_t base =
      std::max(edge_count, g_ADGraph->intermediateEdgeSlotRegistry.size());

  // Reverse propagation creates additional intermediate directional edges
  // beyond the initial soEdges snapshot. Over-allocate to avoid repeated
  // EnsureSlotsPreserve() realloc/copy in the hot loop.
  return base * 4 + 1024;
}

'''
    s = s[:anchor] + helper + s[anchor:]

start, end = find_function(s, "EnableFlatIntermediateDirectionalBackend")
body = s[start:end]

old = '''  g_ADGraph->flatIntermediateDirectionalValues.Resize(
      static_cast<size_t>(g_ADGraph->nBatchDirections),
      g_ADGraph->intermediateEdgeSlotRegistry.size());'''

new = '''  g_ADGraph->flatIntermediateDirectionalValues.Resize(
      static_cast<size_t>(g_ADGraph->nBatchDirections),
      EstimateFlatIntermediateDirectionalSlotCapacity());'''

if old in body:
    body = body.replace(old, new, 1)
    s = s[:start] + body + s[end:]
elif "EstimateFlatIntermediateDirectionalSlotCapacity()" in body:
    print("EnableFlatIntermediateDirectionalBackend already preallocates.")
else:
    raise SystemExit("Could not find flatIntermediateDirectionalValues.Resize block")

p.write_text(s)
PYEOF

python3 /tmp/preallocate_flat_slots.py

cat <<'EOF'

Installed flat intermediate preallocation patch.

Run:
  ./run_had_quadra_nonzero_batch_directional_test.sh
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

Watch:
  - grad diff = 0
  - obj diff = 0
  - reverse time
  - EnsureSlotsPreserve samples should mostly disappear

EOF
