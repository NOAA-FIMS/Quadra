#!/usr/bin/env bash
set -euo pipefail

mkdir -p .quadra_patch_backups

target="core/had_quadra.hpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

if [[ ! -f core/had/batch_directional_flat_accumulator.hpp ]]; then
  echo "ERROR: missing core/had/batch_directional_flat_accumulator.hpp"
  exit 1
fi

if [[ ! -f core/had/intermediate_edge_slot_registry.hpp ]]; then
  echo "ERROR: missing core/had/intermediate_edge_slot_registry.hpp"
  exit 1
fi

cp "$target" ".quadra_patch_backups/had_quadra.hpp.flat_intermediate_backend.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/flat_intermediate_backend.py <<'PYEOF'
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

# Includes
if "batch_directional_flat_accumulator.hpp" not in s:
    s = s.replace("#include <vector>", '#include <vector>\n#include "had/batch_directional_flat_accumulator.hpp"', 1)
if "intermediate_edge_slot_registry.hpp" not in s:
    s = s.replace("#include <vector>", '#include <vector>\n#include "had/intermediate_edge_slot_registry.hpp"', 1)

# Ensure ADGraph registry fields exist.
if "intermediateEdgeSlotRegistry" not in s:
    needle = "std::vector<BTree> soEdgesDot;"
    idx = s.find(needle)
    if idx < 0:
        raise SystemExit("soEdgesDot field not found")
    s = s[:idx+len(needle)] + '''
  had::IntermediateEdgeSlotRegistry intermediateEdgeSlotRegistry;
  bool intermediateEdgeSlotRegistryBuilt = false;
''' + s[idx+len(needle):]

if "flatIntermediateDirectionalValues" not in s:
    needle = "bool intermediateEdgeSlotRegistryBuilt = false;"
    idx = s.find(needle)
    if idx < 0:
        raise SystemExit("intermediateEdgeSlotRegistryBuilt field not found")
    s = s[:idx+len(needle)] + '''
  bool useFlatIntermediateDirectionalBackend = false;
  had::BatchDirectionalFlatAccumulator flatIntermediateDirectionalValues;
''' + s[idx+len(needle):]

# Ensure BuildIntermediateEdgeSlotRegistryFromSoEdges exists.
if "BuildIntermediateEdgeSlotRegistryFromSoEdges" not in s:
    anchor = s.find("inline void ResizeDirectionalBatch")
    helper = '''
struct IntermediateEdgeSlotRegistryDiagnostic {
  std::size_t slots = 0;
  std::size_t source_edges = 0;
};

inline IntermediateEdgeSlotRegistryDiagnostic
BuildIntermediateEdgeSlotRegistryFromSoEdges() {
  IntermediateEdgeSlotRegistryDiagnostic out;
  g_ADGraph->intermediateEdgeSlotRegistry.Clear();

  for (VertexId outer = 0;
       outer < static_cast<VertexId>(g_ADGraph->soEdges.size());
       ++outer) {
    auto &tree = g_ADGraph->soEdges[outer];
    for (const auto &node : tree.nodes) {
      g_ADGraph->intermediateEdgeSlotRegistry.GetOrCreate(
          outer, static_cast<VertexId>(node.key));
      ++out.source_edges;
    }
  }

  g_ADGraph->intermediateEdgeSlotRegistryBuilt = true;
  out.slots = g_ADGraph->intermediateEdgeSlotRegistry.size();
  return out;
}
'''
    if anchor < 0:
        raise SystemExit("ResizeDirectionalBatch not found")
    s = s[:anchor] + helper + s[anchor:]

# Add flat helpers.
if "EnableFlatIntermediateDirectionalBackend" not in s:
    anchor = s.find("inline void ResizeDirectionalBatch")
    helper = '''
inline void ClearFlatIntermediateDirectionalValues() {
  if (g_ADGraph->useFlatIntermediateDirectionalBackend) {
    g_ADGraph->flatIntermediateDirectionalValues.Clear();
  }
}

inline void EnableFlatIntermediateDirectionalBackend() {
  if (!g_ADGraph->intermediateEdgeSlotRegistryBuilt ||
      g_ADGraph->intermediateEdgeSlotRegistry.size() == 0) {
    BuildIntermediateEdgeSlotRegistryFromSoEdges();
  }

  g_ADGraph->useFlatIntermediateDirectionalBackend = true;
  g_ADGraph->flatIntermediateDirectionalValues.Resize(
      static_cast<size_t>(g_ADGraph->nBatchDirections),
      g_ADGraph->intermediateEdgeSlotRegistry.size());
}

inline bool AddFlatIntermediateDirectionalValue(const size_t direction,
                                                const VertexId i,
                                                const VertexId j,
                                                const Real value) {
  if (!g_ADGraph->useFlatIntermediateDirectionalBackend) {
    return false;
  }

  std::size_t slot = 0;
  if (!g_ADGraph->intermediateEdgeSlotRegistry.TryGet(i, j, slot)) {
    return false;
  }

  g_ADGraph->flatIntermediateDirectionalValues.Add(direction, slot, value);
  return true;
}

inline bool TryGetFlatIntermediateDirectionalValue(const size_t direction,
                                                   const VertexId i,
                                                   const VertexId j,
                                                   Real &value_out) {
  if (!g_ADGraph->useFlatIntermediateDirectionalBackend) {
    return false;
  }

  std::size_t slot = 0;
  if (!g_ADGraph->intermediateEdgeSlotRegistry.TryGet(i, j, slot)) {
    return false;
  }

  value_out = g_ADGraph->flatIntermediateDirectionalValues(direction, slot);
  return true;
}
'''
    if anchor < 0:
        raise SystemExit("ResizeDirectionalBatch not found for helpers")
    s = s[:anchor] + helper + s[anchor:]

# Insert clear/enable into PropagateAdjointDirectionalBatch.
start, end = find_function(s, "PropagateAdjointDirectionalBatch")
body = s[start:end]

if "ClearFlatIntermediateDirectionalValues();" not in body:
    if "ClearBatchDirectionalSlotValues();" in body:
        body = body.replace("ClearBatchDirectionalSlotValues();",
                            "ClearBatchDirectionalSlotValues();\n  ClearFlatIntermediateDirectionalValues();",
                            1)
    else:
        brace = body.find("{")
        body = body[:brace+1] + "\n  ClearFlatIntermediateDirectionalValues();\n" + body[brace+1:]

if "EnableFlatIntermediateDirectionalBackend();" not in body:
    marker = "ComputeBatchActiveDirectionMasks(nDirections);"
    if marker in body:
        body = body.replace(marker, "EnableFlatIntermediateDirectionalBackend();\n\n  " + marker, 1)
    else:
        marker = "for (VertexId vid = n_vertices - 1; vid > 0; --vid)"
        if marker in body:
            body = body.replace(marker, "EnableFlatIntermediateDirectionalBackend();\n\n  " + marker, 1)

# Replace query pattern.
old_query = '''        BTree &btreeDot = g_ADGraph->soEdgesDotBatch[static_cast<size_t>(k)][vid];
        ++g_batch_query_count;
        const Real soDot = btreeDot.Query(it->key);'''
new_query = '''        Real soDot = Real(0.0);
        if (!TryGetFlatIntermediateDirectionalValue(
                static_cast<size_t>(k), vid, static_cast<VertexId>(it->key), soDot)) {
          BTree &btreeDot = g_ADGraph->soEdgesDotBatch[static_cast<size_t>(k)][vid];
          ++g_batch_query_count;
          soDot = btreeDot.Query(it->key);
        }'''
if old_query in body:
    body = body.replace(old_query, new_query)
else:
    print("WARNING: did not find primary query pattern")

# Replace cross/create dot insertion patterns.
body = body.replace(
'''            EnsureBatchDotTreeSlot(
                kk, std::max(e1.to, e2.to), "crossDot")
                .Insert(std::min(e1.to, e2.to), crossDot);''',
'''            if (!AddFlatIntermediateDirectionalValue(
                    kk, std::max(e1.to, e2.to), std::min(e1.to, e2.to), crossDot)) {
              EnsureBatchDotTreeSlot(
                  kk, std::max(e1.to, e2.to), "crossDot")
                  .Insert(std::min(e1.to, e2.to), crossDot);
            }'''
)

body = body.replace(
'''            EnsureBatchDotTreeSlot(
                kk, std::max(e1.to, e2.to), "createDot")
                .Insert(std::min(e1.to, e2.to), createDot);''',
'''            if (!AddFlatIntermediateDirectionalValue(
                    kk, std::max(e1.to, e2.to), std::min(e1.to, e2.to), createDot)) {
              EnsureBatchDotTreeSlot(
                  kk, std::max(e1.to, e2.to), "createDot")
                  .Insert(std::min(e1.to, e2.to), createDot);
            }'''
)

s = s[:start] + body + s[end:]

# Patch PushEdgeDotBatchValue off-diagonal fallback insert.
start, end = find_function(s, "PushEdgeDotBatchValue")
body = s[start:end]
if "AddFlatIntermediateDirectionalValue" not in body:
    old = "      trees[outer].Insert(inner, valDot);"
    new = '''      if (!AddFlatIntermediateDirectionalValue(direction, outer, inner, valDot)) {
        trees[outer].Insert(inner, valDot);
      }'''
    if old in body:
        body = body.replace(old, new, 1)
    else:
        print("WARNING: did not find PushEdgeDotBatchValue insert pattern")
    s = s[:start] + body + s[end:]

p.write_text(s)
PYEOF

python3 /tmp/flat_intermediate_backend.py

cat <<'EOF'

Installed experimental flat intermediate directional backend.

Run:
  ./run_had_quadra_nonzero_batch_directional_test.sh
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

EOF
