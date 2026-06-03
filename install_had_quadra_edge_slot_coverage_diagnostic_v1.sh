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

cp "$header" ".quadra_patch_backups/had_quadra.hpp.edge_slot_coverage_diag.$(date +%Y%m%d_%H%M%S).bak"
cp "$bench" ".quadra_patch_backups/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp.edge_slot_coverage_diag.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/edge_slot_coverage_diag.py <<'PYEOF'
from pathlib import Path

h = Path("core/had_quadra.hpp")
s = h.read_text()

if "CountMappedKeys" not in s:
    idx = s.find("inline Real Query")
    if idx < 0:
        raise SystemExit("BTree Query not found")

    brace = s.find("{", idx)
    depth = 0
    end = None
    for i in range(brace, len(s)):
        if s[i] == "{":
            depth += 1
        elif s[i] == "}":
            depth -= 1
            if depth == 0:
                end = i + 1
                break
    if end is None:
        raise SystemExit("Could not find end of Query")

    helper = """

  inline std::size_t CountMappedKeys(const BTree &slot_map_tree) const {
    std::size_t mapped = 0;
    for (const auto &node : nodes) {
      if (slot_map_tree.Query(node.key) != Real(0.0)) {
        ++mapped;
      }
    }
    return mapped;
  }
"""
    s = s[:end] + helper + s[end:]

if "BatchEdgeSlotCoverageDiagnostic" not in s:
    anchor = s.find("inline void ResizeDirectionalBatch")
    if anchor < 0:
        raise SystemExit("ResizeDirectionalBatch not found")

    diag = """
struct BatchEdgeSlotCoverageDiagnostic {
  std::size_t total_edges = 0;
  std::size_t mapped_edges = 0;
  std::size_t unmapped_edges = 0;
  double mapped_fraction = 0.0;
};

inline BatchEdgeSlotCoverageDiagnostic DiagnoseBatchEdgeSlotCoverage() {
  BatchEdgeSlotCoverageDiagnostic out;

  if (!g_ADGraph->useBatchDirectionalSlotWorkspace) {
    return out;
  }

  const std::size_t n =
      std::min(g_ADGraph->soEdges.size(),
               g_ADGraph->batchSlotOuterInnerToSlot.size());

  for (std::size_t outer = 0; outer < n; ++outer) {
    const auto &edges = g_ADGraph->soEdges[outer];
    const auto &slot_map = g_ADGraph->batchSlotOuterInnerToSlot[outer];

    const std::size_t total = edges.Size();
    const std::size_t mapped = edges.CountMappedKeys(slot_map);

    out.total_edges += total;
    out.mapped_edges += mapped;
  }

  out.unmapped_edges = out.total_edges - out.mapped_edges;
  if (out.total_edges > 0) {
    out.mapped_fraction =
        static_cast<double>(out.mapped_edges) /
        static_cast<double>(out.total_edges);
  }

  return out;
}

inline void PrintBatchEdgeSlotCoverageDiagnostic() {
  const auto d = DiagnoseBatchEdgeSlotCoverage();
  std::cerr << "batch edge slot coverage:"
            << " total=" << d.total_edges
            << " mapped=" << d.mapped_edges
            << " unmapped=" << d.unmapped_edges
            << " mapped_fraction=" << d.mapped_fraction
            << "\\n";
}

"""
    s = s[:anchor] + diag + s[anchor:]

h.write_text(s)

b = Path("benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp")
t = b.read_text()

old = """        workspace.propagate_directional_batch();

        const auto rb3 = Clock::now();"""

new = """        workspace.propagate_directional_batch();

        if (m == 500 && r == reps - 1) {
            had::PrintBatchEdgeSlotCoverageDiagnostic();
        }

        const auto rb3 = Clock::now();"""

if old not in t:
    raise SystemExit("Could not find benchmark directional propagation block")

t = t.replace(old, new, 1)
b.write_text(t)
PYEOF

python3 /tmp/edge_slot_coverage_diag.py

cat <<'EOF'

Installed edge-slot coverage diagnostic.

Run:
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

Look for:
  batch edge slot coverage: total=... mapped=... unmapped=... mapped_fraction=...

EOF
