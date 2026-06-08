#!/usr/bin/env bash
set -euo pipefail

# install_flat_intermediate_missing_key_samples_v1.sh
#
# Diagnostics for zero flat-intermediate hits:
#
#   flat intermediate: read_hit=0 read_miss=42495 write_hit=0 write_miss=59990
#
# This means the flat backend is enabled, but the registry does not contain
# the (outer, inner) keys being queried/written.
#
# This patch records a few missing read/write keys so we can see whether:
#   - registry is built too early/late
#   - keys are using intermediate AD vertex ids not present in soEdges snapshot
#   - write/read arguments are reversed or otherwise mismatched
#
# Diagnostic-only.

mkdir -p .quadra_patch_backups

target="core/had_quadra.hpp"
bench="benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi
if [[ ! -f "$bench" ]]; then
  echo "ERROR: missing $bench"
  exit 1
fi

cp "$target" ".quadra_patch_backups/had_quadra.hpp.flat_missing_keys.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/flat_missing_keys.py <<'PYEOF'
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

# Add sample globals/helpers after flat counters.
if "g_flat_intermediate_read_miss_samples" not in s:
    anchor = "inline std::uint64_t g_flat_intermediate_write_miss_count = 0;"
    if anchor not in s:
        raise SystemExit("flat intermediate counters not found")

    insertion = """

inline std::vector<std::pair<VertexId, VertexId>>
    g_flat_intermediate_read_miss_samples;
inline std::vector<std::pair<VertexId, VertexId>>
    g_flat_intermediate_write_miss_samples;

inline void RecordFlatIntermediateMissSample(
    std::vector<std::pair<VertexId, VertexId>> &samples,
    const VertexId i,
    const VertexId j) {
  if (samples.size() < 12) {
    samples.emplace_back(std::max(i, j), std::min(i, j));
  }
}

inline void PrintFlatIntermediateMissSamples() {
  std::cerr << "flat intermediate read miss samples:";
  for (const auto &p : g_flat_intermediate_read_miss_samples) {
    std::cerr << " (" << p.first << "," << p.second << ")";
  }
  std::cerr << "\\n";

  std::cerr << "flat intermediate write miss samples:";
  for (const auto &p : g_flat_intermediate_write_miss_samples) {
    std::cerr << " (" << p.first << "," << p.second << ")";
  }
  std::cerr << "\\n";
}
"""
    s = s.replace(anchor, anchor + insertion, 1)

# Reset sample vectors in ResetFlatIntermediateDirectionalCounters.
start, end = find_function(s, "ResetFlatIntermediateDirectionalCounters")
body = s[start:end]
if "g_flat_intermediate_read_miss_samples.clear();" not in body:
    close = body.rfind("}")
    body = body[:close] + """
  g_flat_intermediate_read_miss_samples.clear();
  g_flat_intermediate_write_miss_samples.clear();
""" + body[close:]
    s = s[:start] + body + s[end:]

# Print samples in PrintFlatIntermediateDirectionalCounters.
start, end = find_function(s, "PrintFlatIntermediateDirectionalCounters")
body = s[start:end]
if "PrintFlatIntermediateMissSamples();" not in body:
    close = body.rfind("}")
    body = body[:close] + """
  PrintFlatIntermediateMissSamples();
""" + body[close:]
    s = s[:start] + body + s[end:]

# Record write misses.
start, end = find_function(s, "AddFlatIntermediateDirectionalValue")
body = s[start:end]
if "g_flat_intermediate_write_miss_samples" not in body:
    body = body.replace(
        """  if (!g_ADGraph->useFlatIntermediateDirectionalBackend) {
    ++g_flat_intermediate_write_miss_count;
    return false;
  }""",
        """  if (!g_ADGraph->useFlatIntermediateDirectionalBackend) {
    ++g_flat_intermediate_write_miss_count;
    RecordFlatIntermediateMissSample(
        g_flat_intermediate_write_miss_samples, i, j);
    return false;
  }""",
        1,
    )

    body = body.replace(
        """  if (!g_ADGraph->intermediateEdgeSlotRegistry.TryGet(i, j, slot)) {
    ++g_flat_intermediate_write_miss_count;
    return false;
  }""",
        """  if (!g_ADGraph->intermediateEdgeSlotRegistry.TryGet(i, j, slot)) {
    ++g_flat_intermediate_write_miss_count;
    RecordFlatIntermediateMissSample(
        g_flat_intermediate_write_miss_samples, i, j);
    return false;
  }""",
        1,
    )

s = s[:start] + body + s[end:]

# Record read misses.
start, end = find_function(s, "TryGetFlatIntermediateDirectionalValue")
body = s[start:end]
if "g_flat_intermediate_read_miss_samples" not in body:
    body = body.replace(
        """  if (!g_ADGraph->useFlatIntermediateDirectionalBackend) {
    ++g_flat_intermediate_read_miss_count;
    return false;
  }""",
        """  if (!g_ADGraph->useFlatIntermediateDirectionalBackend) {
    ++g_flat_intermediate_read_miss_count;
    RecordFlatIntermediateMissSample(
        g_flat_intermediate_read_miss_samples, i, j);
    return false;
  }""",
        1,
    )

    body = body.replace(
        """  if (!g_ADGraph->intermediateEdgeSlotRegistry.TryGet(i, j, slot)) {
    ++g_flat_intermediate_read_miss_count;
    return false;
  }""",
        """  if (!g_ADGraph->intermediateEdgeSlotRegistry.TryGet(i, j, slot)) {
    ++g_flat_intermediate_read_miss_count;
    RecordFlatIntermediateMissSample(
        g_flat_intermediate_read_miss_samples, i, j);
    return false;
  }""",
        1,
    )

s = s[:start] + body + s[end:]

p.write_text(s)
PYEOF

python3 /tmp/flat_missing_keys.py

cat <<'EOF'

Installed flat intermediate missing-key samples.

Run:
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

Look for:
  flat intermediate read miss samples: ...
  flat intermediate write miss samples: ...

EOF
