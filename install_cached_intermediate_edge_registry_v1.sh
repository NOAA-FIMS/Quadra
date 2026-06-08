#!/usr/bin/env bash
set -euo pipefail

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

if [[ ! -f core/had/intermediate_edge_slot_registry.hpp ]]; then
  echo "ERROR: missing core/had/intermediate_edge_slot_registry.hpp"
  echo "Run install_intermediate_edge_slot_registry_scaffold_v1.sh first."
  exit 1
fi

cp "$header" ".quadra_patch_backups/had_quadra.hpp.cached_intermediate_registry.$(date +%Y%m%d_%H%M%S).bak"
cp "$bench" ".quadra_patch_backups/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp.cached_intermediate_registry.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/cached_intermediate_registry.py <<'PYEOF'
from pathlib import Path

h = Path("core/had_quadra.hpp")
s = h.read_text()

if "intermediate_edge_slot_registry.hpp" not in s:
    if "#include <vector>" in s:
        s = s.replace(
            "#include <vector>",
            '#include <vector>\n#include "had/intermediate_edge_slot_registry.hpp"',
            1,
        )
    else:
        s = '#include "had/intermediate_edge_slot_registry.hpp"\n' + s

if "intermediateEdgeSlotRegistry" not in s:
    needle = "std::vector<BTree> soEdgesDot;"
    idx = s.find(needle)
    if idx < 0:
        raise SystemExit("Could not find soEdgesDot field")
    insert_at = idx + len(needle)
    insertion = """
  // Cached slot registry for intermediate second-order edges.
  // Future flat reverse backend will use this to replace BTree query/insert
  // with direct slot-indexed directional accumulation.
  had::IntermediateEdgeSlotRegistry intermediateEdgeSlotRegistry;
  bool intermediateEdgeSlotRegistryBuilt = false;"""
    s = s[:insert_at] + insertion + s[insert_at:]

if "intermediateEdgeSlotRegistry.Clear();" not in s:
    needle = "soEdgesDot.clear();"
    if needle in s:
        s = s.replace(
            needle,
            needle + """
    intermediateEdgeSlotRegistry.Clear();
    intermediateEdgeSlotRegistryBuilt = false;""",
            1,
        )

if "BuildIntermediateEdgeSlotRegistryFromSoEdges" not in s:
    anchor = s.find("inline void ResizeDirectionalBatch")
    if anchor < 0:
        raise SystemExit("ResizeDirectionalBatch not found")

    helper = """
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
          outer,
          static_cast<VertexId>(node.key));
      ++out.source_edges;
    }
  }

  g_ADGraph->intermediateEdgeSlotRegistryBuilt = true;
  out.slots = g_ADGraph->intermediateEdgeSlotRegistry.size();
  return out;
}

inline IntermediateEdgeSlotRegistryDiagnostic
GetIntermediateEdgeSlotRegistryDiagnostic() {
  IntermediateEdgeSlotRegistryDiagnostic out;
  out.slots = g_ADGraph->intermediateEdgeSlotRegistry.size();
  out.source_edges = g_ADGraph->intermediateEdgeSlotRegistry.size();
  return out;
}

inline void PrintIntermediateEdgeSlotRegistryDiagnostic() {
  const auto d = BuildIntermediateEdgeSlotRegistryFromSoEdges();
  std::cerr << "intermediate edge registry:"
            << " source_edges=" << d.source_edges
            << " slots=" << d.slots
            << "\\n";
}

"""
    s = s[:anchor] + helper + s[anchor:]

h.write_text(s)

b = Path("benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp")
t = b.read_text()

old = """        workspace.propagate_directional_batch();

        const auto rb3 = Clock::now();"""

new = """        workspace.propagate_directional_batch();

        if (m == 500 && r == reps - 1) {
            had::PrintIntermediateEdgeSlotRegistryDiagnostic();
        }

        const auto rb3 = Clock::now();"""

if old not in t:
    raise SystemExit("Could not find propagation timing block in benchmark")

t = t.replace(old, new, 1)
b.write_text(t)
PYEOF

python3 /tmp/cached_intermediate_registry.py

cat <<'EOF'

Installed cached intermediate edge registry scaffold.

Run:
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

Look for:
  intermediate edge registry: source_edges=... slots=...

EOF
