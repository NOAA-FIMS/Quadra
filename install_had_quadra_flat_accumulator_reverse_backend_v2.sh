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
  echo "Run install_batch_directional_flat_accumulator_scaffold_v1.sh first."
  exit 1
fi

cp "$target" ".quadra_patch_backups/had_quadra.hpp.flat_accumulator_reverse_backend_v2.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/flat_accumulator_reverse_backend_v2.py <<'PYEOF'
from pathlib import Path

p = Path("core/had_quadra.hpp")
s = p.read_text()

def find_function(src, name):
    idx = src.find(name)
    if idx < 0:
        raise SystemExit(f"function {name} not found")
    start = src.rfind("\n", 0, idx) + 1
    brace = src.find("{", idx)
    if brace < 0:
        raise SystemExit(f"function {name} has no opening brace")
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
        raise SystemExit(f"function {name} has no closing brace")
    return start, end

if "batch_directional_flat_accumulator.hpp" not in s:
    if "#include <vector>" in s:
        s = s.replace("#include <vector>",
                      '#include <vector>\n#include "had/batch_directional_flat_accumulator.hpp"',
                      1)
    else:
        s = '#include "had/batch_directional_flat_accumulator.hpp"\n' + s

if "useBatchDirectionalFlatAccumulator" not in s:
    needle = "std::vector<std::vector<Real>> batchDirectionalSlotValues;"
    idx = s.find(needle)
    if idx < 0:
        raise SystemExit("batchDirectionalSlotValues field not found")
    insert_at = idx + len(needle)
    s = s[:insert_at] + """

        // Optional flat accumulator backend for mapped directional Hdot slots.
        bool useBatchDirectionalFlatAccumulator = false;
        had::BatchDirectionalFlatAccumulator batchDirectionalFlatAccumulator;""" + s[insert_at:]

if "EnableBatchDirectionalFlatAccumulator" not in s:
    anchor = s.find("inline void DisableBatchDirectionalSlotWorkspace()")
    if anchor < 0:
        raise SystemExit("DisableBatchDirectionalSlotWorkspace not found")
    helper = """
inline void EnableBatchDirectionalFlatAccumulator() {
  if (!g_ADGraph->useBatchDirectionalSlotWorkspace) {
    throw std::runtime_error(
        "EnableBatchDirectionalFlatAccumulator requires slot workspace enabled");
  }

  const size_t nDirections =
      static_cast<size_t>(g_ADGraph->nBatchDirections);

  const size_t nSlots =
      g_ADGraph->batchDirectionalSlotValues.empty()
          ? 0
          : g_ADGraph->batchDirectionalSlotValues.front().size();

  g_ADGraph->useBatchDirectionalFlatAccumulator = true;
  g_ADGraph->batchDirectionalFlatAccumulator.Resize(nDirections, nSlots);
}

inline void DisableBatchDirectionalFlatAccumulator() {
  g_ADGraph->useBatchDirectionalFlatAccumulator = false;
  g_ADGraph->batchDirectionalFlatAccumulator =
      had::BatchDirectionalFlatAccumulator();
}

"""
    s = s[:anchor] + helper + s[anchor:]

start, end = find_function(s, "DisableBatchDirectionalSlotWorkspace")
body = s[start:end]
if "useBatchDirectionalFlatAccumulator" not in body:
    brace = body.find("{")
    body = body[:brace+1] + """
  g_ADGraph->useBatchDirectionalFlatAccumulator = false;
  g_ADGraph->batchDirectionalFlatAccumulator =
      had::BatchDirectionalFlatAccumulator();
""" + body[brace+1:]
    s = s[:start] + body + s[end:]

new_add = """
inline bool AddBatchDirectionalSlotValue(const size_t direction,
                                         const VertexId i,
                                         const VertexId j,
                                         const Real value) {
  if (!g_ADGraph->useBatchDirectionalSlotWorkspace) {
    ++g_batch_slot_write_miss_count;
    return false;
  }
  if (direction >= g_ADGraph->batchDirectionalSlotValues.size()) {
    ++g_batch_slot_write_miss_count;
    return false;
  }

  int slot = -1;

  if (i == j) {
    if (i < g_ADGraph->batchSelfSlot.size()) {
      slot = g_ADGraph->batchSelfSlot[i];
    }
  } else {
    const VertexId outer = std::max(i, j);
    const VertexId inner = std::min(i, j);

    if (outer < g_ADGraph->batchSlotOuterInnerToSlot.size()) {
      const Real stored =
          g_ADGraph->batchSlotOuterInnerToSlot[outer].Query(inner);
      if (stored != Real(0.0)) {
        slot = static_cast<int>(stored) - 1;
      }
    }
  }

  if (slot < 0) {
    ++g_batch_slot_write_miss_count;
    return false;
  }

  auto &values = g_ADGraph->batchDirectionalSlotValues[direction];
  if (slot >= static_cast<int>(values.size())) {
    ++g_batch_slot_write_miss_count;
    return false;
  }

  if (g_ADGraph->useBatchDirectionalFlatAccumulator) {
    g_ADGraph->batchDirectionalFlatAccumulator.Add(
        direction, static_cast<size_t>(slot), value);
  } else {
    values[static_cast<size_t>(slot)] += value;
  }

  ++g_batch_slot_write_hit_count;
  return true;
}
"""
start, end = find_function(s, "AddBatchDirectionalSlotValue")
s = s[:start] + new_add + s[end:]

new_get = """
inline Real GetBatchDirectionalSlotValue(const size_t direction,
                                         const int slot) {
  if (g_ADGraph->useBatchDirectionalFlatAccumulator) {
    return g_ADGraph->batchDirectionalFlatAccumulator(
        direction, static_cast<size_t>(slot));
  }

  return g_ADGraph->batchDirectionalSlotValues[direction][static_cast<size_t>(slot)];
}
"""
start, end = find_function(s, "GetBatchDirectionalSlotValue")
s = s[:start] + new_get + s[end:]

new_try = """
inline bool TryGetBatchDirectionalSlotValue(const size_t direction,
                                            const VertexId i,
                                            const VertexId j,
                                            Real &value_out) {
  int slot = -1;
  if (!TryGetBatchDirectionalSlot(i, j, slot)) {
    ++g_batch_slot_query_miss_count;
    return false;
  }

  if (direction >= g_ADGraph->batchDirectionalSlotValues.size()) {
    ++g_batch_slot_query_miss_count;
    return false;
  }

  const auto &values = g_ADGraph->batchDirectionalSlotValues[direction];

  if (slot < 0 || slot >= static_cast<int>(values.size())) {
    ++g_batch_slot_query_miss_count;
    return false;
  }

  if (g_ADGraph->useBatchDirectionalFlatAccumulator) {
    value_out = g_ADGraph->batchDirectionalFlatAccumulator(
        direction, static_cast<size_t>(slot));
  } else {
    value_out = values[static_cast<size_t>(slot)];
  }

  ++g_batch_slot_query_hit_count;
  return true;
}
"""
start, end = find_function(s, "TryGetBatchDirectionalSlotValue")
s = s[:start] + new_try + s[end:]

start, end = find_function(s, "EnableBatchDirectionalSlotWorkspace")
body = s[start:end]
if "batchDirectionalFlatAccumulator.Resize" not in body:
    close = body.rfind("}")
    body = body[:close] + """
  g_ADGraph->useBatchDirectionalFlatAccumulator = true;
  g_ADGraph->batchDirectionalFlatAccumulator.Resize(
      static_cast<size_t>(nDirections),
      static_cast<size_t>(nSlots));
""" + body[close:]
    s = s[:start] + body + s[end:]

start, end = find_function(s, "ClearBatchDirectionalSlotValues")
body = s[start:end]
if "batchDirectionalFlatAccumulator.Clear" not in body:
    brace = body.find("{")
    body = body[:brace+1] + """
  if (g_ADGraph->useBatchDirectionalFlatAccumulator) {
    g_ADGraph->batchDirectionalFlatAccumulator.Clear();
  }
""" + body[brace+1:]
    s = s[:start] + body + s[end:]

start, end = find_function(s, "ResizeDirectionalBatch")
body = s[start:end]
if "batchDirectionalFlatAccumulator.Resize" not in body:
    close = body.rfind("}")
    body = body[:close] + """
        if (g_ADGraph->useBatchDirectionalFlatAccumulator)
        {
            const size_t nSlots =
                g_ADGraph->batchDirectionalSlotValues.empty()
                    ? 0
                    : g_ADGraph->batchDirectionalSlotValues.front().size();
            g_ADGraph->batchDirectionalFlatAccumulator.Resize(
                static_cast<size_t>(nDirections), nSlots);
        }
""" + body[close:]
    s = s[:start] + body + s[end:]

p.write_text(s)
PYEOF

python3 /tmp/flat_accumulator_reverse_backend_v2.py

cat <<'EOF'

Installed experimental flat accumulator reverse backend v2.

Run:
  ./run_had_quadra_nonzero_batch_directional_test.sh
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10
  ./run_sparse_rw1_flat_accumulator_trace_check.sh 10

EOF
