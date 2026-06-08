#!/usr/bin/env bash
set -euo pipefail

mkdir -p .quadra_patch_backups

target="core/had_quadra.hpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/had_quadra.hpp.flat_accumulator_reverse_backend.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/flat_accumulator_reverse_backend.py <<'PYEOF'
from pathlib import Path

p = Path("core/had_quadra.hpp")
s = p.read_text()

if "batch_directional_flat_accumulator.hpp" not in s:
    if "#include <vector>" in s:
        s = s.replace("#include <vector>", '#include <vector>\n#include "had/batch_directional_flat_accumulator.hpp"', 1)
    else:
        s = '#include "had/batch_directional_flat_accumulator.hpp"\n' + s

if "useBatchDirectionalFlatAccumulator" not in s:
    old = "        std::vector<std::vector<Real>> batchDirectionalSlotValues;"
    new = '''        std::vector<std::vector<Real>> batchDirectionalSlotValues;

        // Optional flat accumulator backend for mapped directional Hdot slots.
        bool useBatchDirectionalFlatAccumulator = false;
        had::BatchDirectionalFlatAccumulator batchDirectionalFlatAccumulator;'''
    if old not in s:
        raise SystemExit("Could not find batchDirectionalSlotValues field")
    s = s.replace(old, new, 1)

if "batchDirectionalFlatAccumulator = had::BatchDirectionalFlatAccumulator();" not in s:
    old = "    batchDirectionalSlotValues.clear();"
    new = '''    batchDirectionalSlotValues.clear();
    useBatchDirectionalFlatAccumulator = false;
    batchDirectionalFlatAccumulator = had::BatchDirectionalFlatAccumulator();'''
    if old in s:
        s = s.replace(old, new, 1)

if "EnableBatchDirectionalFlatAccumulator" not in s:
    anchor = s.find("inline void DisableBatchDirectionalSlotWorkspace()")
    if anchor < 0:
        raise SystemExit("DisableBatchDirectionalSlotWorkspace not found")

    helper = '''
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

'''
    s = s[:anchor] + helper + s[anchor:]

old = '''inline void DisableBatchDirectionalSlotWorkspace() {
  g_ADGraph->useBatchDirectionalSlotWorkspace = false;'''
new = '''inline void DisableBatchDirectionalSlotWorkspace() {
  g_ADGraph->useBatchDirectionalFlatAccumulator = false;
  g_ADGraph->batchDirectionalFlatAccumulator =
      had::BatchDirectionalFlatAccumulator();
  g_ADGraph->useBatchDirectionalSlotWorkspace = false;'''
s = s.replace(old, new, 1)

old = '''  g_ADGraph->batchDirectionalSlotValues.assign(
      static_cast<size_t>(nDirections),
      std::vector<Real>(static_cast<size_t>(nSlots), Real(0.0)));
}'''
new = '''  g_ADGraph->batchDirectionalSlotValues.assign(
      static_cast<size_t>(nDirections),
      std::vector<Real>(static_cast<size_t>(nSlots), Real(0.0)));

  g_ADGraph->useBatchDirectionalFlatAccumulator = true;
  g_ADGraph->batchDirectionalFlatAccumulator.Resize(
      static_cast<size_t>(nDirections),
      static_cast<size_t>(nSlots));
}'''
if old in s:
    s = s.replace(old, new, 1)
else:
    print("WARNING: EnableBatchDirectionalSlotWorkspace assignment block not found")

old = '''  for (auto &values : g_ADGraph->batchDirectionalSlotValues) {
    std::fill(values.begin(), values.end(), Real(0.0));
  }
}'''
new = '''  if (g_ADGraph->useBatchDirectionalFlatAccumulator) {
    g_ADGraph->batchDirectionalFlatAccumulator.Clear();
  }

  for (auto &values : g_ADGraph->batchDirectionalSlotValues) {
    std::fill(values.begin(), values.end(), Real(0.0));
  }
}'''
s = s.replace(old, new, 1)

old = '''            g_ADGraph->batchDirectionalSlotValues.assign(
                static_cast<size_t>(nDirections),
                std::vector<Real>(nSlots, Real(0.0)));
        }'''
new = '''            g_ADGraph->batchDirectionalSlotValues.assign(
                static_cast<size_t>(nDirections),
                std::vector<Real>(nSlots, Real(0.0)));
            if (g_ADGraph->useBatchDirectionalFlatAccumulator)
            {
                g_ADGraph->batchDirectionalFlatAccumulator.Resize(
                    static_cast<size_t>(nDirections), nSlots);
            }
        }'''
s = s.replace(old, new, 1)

old = '''  values[static_cast<size_t>(slot)] += value;
  ++g_batch_slot_write_hit_count;
  return true;'''
new = '''  if (g_ADGraph->useBatchDirectionalFlatAccumulator) {
    g_ADGraph->batchDirectionalFlatAccumulator.Add(
        direction, static_cast<size_t>(slot), value);
  } else {
    values[static_cast<size_t>(slot)] += value;
  }

  ++g_batch_slot_write_hit_count;
  return true;'''
if old not in s:
    raise SystemExit("Could not find slot write block")
s = s.replace(old, new, 1)

old = '''inline Real GetBatchDirectionalSlotValue(const size_t direction,
                                         const int slot) {
  return g_ADGraph->batchDirectionalSlotValues[direction][static_cast<size_t>(slot)];
}'''
new = '''inline Real GetBatchDirectionalSlotValue(const size_t direction,
                                         const int slot) {
  if (g_ADGraph->useBatchDirectionalFlatAccumulator) {
    return g_ADGraph->batchDirectionalFlatAccumulator(
        direction, static_cast<size_t>(slot));
  }

  return g_ADGraph->batchDirectionalSlotValues[direction][static_cast<size_t>(slot)];
}'''
if old in s:
    s = s.replace(old, new, 1)
else:
    print("WARNING: GetBatchDirectionalSlotValue exact block not found")

old = '''  value_out = values[static_cast<size_t>(slot)];
  ++g_batch_slot_query_hit_count;
  return true;'''
new = '''  if (g_ADGraph->useBatchDirectionalFlatAccumulator) {
    value_out = g_ADGraph->batchDirectionalFlatAccumulator(
        direction, static_cast<size_t>(slot));
  } else {
    value_out = values[static_cast<size_t>(slot)];
  }

  ++g_batch_slot_query_hit_count;
  return true;'''
if old in s:
    s = s.replace(old, new, 1)
else:
    print("WARNING: TryGetBatchDirectionalSlotValue read block not found")

p.write_text(s)
PYEOF

python3 /tmp/flat_accumulator_reverse_backend.py

cat <<'EOF'

Installed experimental flat accumulator reverse backend.

Run:
  ./run_had_quadra_nonzero_batch_directional_test.sh
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10
  ./run_sparse_rw1_flat_accumulator_trace_check.sh 10

EOF
