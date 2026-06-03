#!/usr/bin/env bash
set -euo pipefail

# install_flat_intermediate_backend_hit_diagnostics_v1.sh
#
# Adds diagnostics for the experimental flat intermediate backend:
#   - flat intermediate read hits/misses
#   - flat intermediate write hits/misses
#
# Purpose:
#   The previous patch preserved correctness but did not reduce BTree query count
#   and made reverse slower. We need to determine whether:
#     1. registry slots are not available at read/write time,
#     2. read/write key pairs do not match the registry,
#     3. write path was not patched,
#     4. query count is misleading because fallback still dominates.
#
# This patch is diagnostic-only.

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

cp "$target" ".quadra_patch_backups/had_quadra.hpp.flat_intermediate_hit_diag.$(date +%Y%m%d_%H%M%S).bak"
cp "$bench" ".quadra_patch_backups/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp.flat_intermediate_hit_diag.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/flat_intermediate_hit_diag.py <<'PYEOF'
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

# Add counters before flat helper functions.
if "g_flat_intermediate_read_hit_count" not in s:
    anchor = s.find("inline void ClearFlatIntermediateDirectionalValues")
    if anchor < 0:
        raise SystemExit("ClearFlatIntermediateDirectionalValues not found; flat backend not installed?")

    counters = """
inline std::uint64_t g_flat_intermediate_read_hit_count = 0;
inline std::uint64_t g_flat_intermediate_read_miss_count = 0;
inline std::uint64_t g_flat_intermediate_write_hit_count = 0;
inline std::uint64_t g_flat_intermediate_write_miss_count = 0;

inline void ResetFlatIntermediateDirectionalCounters() {
  g_flat_intermediate_read_hit_count = 0;
  g_flat_intermediate_read_miss_count = 0;
  g_flat_intermediate_write_hit_count = 0;
  g_flat_intermediate_write_miss_count = 0;
}

inline void PrintFlatIntermediateDirectionalCounters() {
  std::cerr << "flat intermediate:"
            << " read_hit=" << g_flat_intermediate_read_hit_count
            << " read_miss=" << g_flat_intermediate_read_miss_count
            << " write_hit=" << g_flat_intermediate_write_hit_count
            << " write_miss=" << g_flat_intermediate_write_miss_count
            << "\\n";
}

"""
    s = s[:anchor] + counters + s[anchor:]

# Reset counters at beginning of PropagateAdjointDirectionalBatch.
start, end = find_function(s, "PropagateAdjointDirectionalBatch")
body = s[start:end]
if "ResetFlatIntermediateDirectionalCounters();" not in body:
    brace = body.find("{")
    body = body[:brace+1] + "\n  ResetFlatIntermediateDirectionalCounters();\n" + body[brace+1:]
    s = s[:start] + body + s[end:]

# Instrument AddFlatIntermediateDirectionalValue.
start, end = find_function(s, "AddFlatIntermediateDirectionalValue")
body = s[start:end]

if "g_flat_intermediate_write_hit_count" not in body:
    body = body.replace(
        """  if (!g_ADGraph->useFlatIntermediateDirectionalBackend) {
    return false;
  }""",
        """  if (!g_ADGraph->useFlatIntermediateDirectionalBackend) {
    ++g_flat_intermediate_write_miss_count;
    return false;
  }""",
        1,
    )

    body = body.replace(
        """  if (!g_ADGraph->intermediateEdgeSlotRegistry.TryGet(i, j, slot)) {
    return false;
  }""",
        """  if (!g_ADGraph->intermediateEdgeSlotRegistry.TryGet(i, j, slot)) {
    ++g_flat_intermediate_write_miss_count;
    return false;
  }""",
        1,
    )

    body = body.replace(
        """  g_ADGraph->flatIntermediateDirectionalValues.Add(direction, slot, value);
  return true;""",
        """  g_ADGraph->flatIntermediateDirectionalValues.Add(direction, slot, value);
  ++g_flat_intermediate_write_hit_count;
  return true;""",
        1,
    )

s = s[:start] + body + s[end:]

# Instrument TryGetFlatIntermediateDirectionalValue.
start, end = find_function(s, "TryGetFlatIntermediateDirectionalValue")
body = s[start:end]

if "g_flat_intermediate_read_hit_count" not in body:
    body = body.replace(
        """  if (!g_ADGraph->useFlatIntermediateDirectionalBackend) {
    return false;
  }""",
        """  if (!g_ADGraph->useFlatIntermediateDirectionalBackend) {
    ++g_flat_intermediate_read_miss_count;
    return false;
  }""",
        1,
    )

    body = body.replace(
        """  if (!g_ADGraph->intermediateEdgeSlotRegistry.TryGet(i, j, slot)) {
    return false;
  }""",
        """  if (!g_ADGraph->intermediateEdgeSlotRegistry.TryGet(i, j, slot)) {
    ++g_flat_intermediate_read_miss_count;
    return false;
  }""",
        1,
    )

    body = body.replace(
        """  value_out = g_ADGraph->flatIntermediateDirectionalValues(direction, slot);
  return true;""",
        """  value_out = g_ADGraph->flatIntermediateDirectionalValues(direction, slot);
  ++g_flat_intermediate_read_hit_count;
  return true;""",
        1,
    )

s = s[:start] + body + s[end:]

p.write_text(s)

# Add one print in benchmark after propagation for m=500 last rep.
b = Path("benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp")
t = b.read_text()

block = """        if (m == 500 && r == reps - 1) {
            had::PrintFlatIntermediateDirectionalCounters();
        }

"""
if "PrintFlatIntermediateDirectionalCounters" not in t:
    marker = """        workspace.propagate_directional_batch();

        const auto rb3 = Clock::now();"""
    replacement = """        workspace.propagate_directional_batch();

        if (m == 500 && r == reps - 1) {
            had::PrintFlatIntermediateDirectionalCounters();
        }

        const auto rb3 = Clock::now();"""
    if marker not in t:
        raise SystemExit("benchmark propagation marker not found")
    t = t.replace(marker, replacement, 1)

b.write_text(t)
PYEOF

python3 /tmp/flat_intermediate_hit_diag.py

cat <<'EOF'

Installed flat intermediate hit/miss diagnostics.

Run:
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

Look for:
  flat intermediate: read_hit=... read_miss=... write_hit=... write_miss=...

EOF
