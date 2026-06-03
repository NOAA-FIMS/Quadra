#!/usr/bin/env bash
set -euo pipefail

mkdir -p .quadra_patch_backups

accum="core/had/batch_directional_flat_accumulator.hpp"
header="core/had_quadra.hpp"

if [[ ! -f "$accum" ]]; then
  echo "ERROR: missing $accum"
  exit 1
fi
if [[ ! -f "$header" ]]; then
  echo "ERROR: missing $header"
  exit 1
fi

cp "$accum" ".quadra_patch_backups/batch_directional_flat_accumulator.hpp.dynamic_slots.$(date +%Y%m%d_%H%M%S).bak"
cp "$header" ".quadra_patch_backups/had_quadra.hpp.dynamic_intermediate_slots.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/repair_dynamic_slots.py <<'PYEOF'
from pathlib import Path

p = Path("core/had/batch_directional_flat_accumulator.hpp")
s = p.read_text()

if "EnsureSlotsPreserve" not in s:
    marker = '''  void Clear() {
    std::fill(values_.begin(), values_.end(), 0.0);
  }

'''
    method = '''  void EnsureSlotsPreserve(std::size_t n_slots) {
    if (n_slots <= n_slots_) {
      return;
    }

    std::vector<double> next(n_directions_ * n_slots, 0.0);

    for (std::size_t k = 0; k < n_directions_; ++k) {
      for (std::size_t slot = 0; slot < n_slots_; ++slot) {
        next[k * n_slots + slot] = values_[k * n_slots_ + slot];
      }
    }

    values_.swap(next);
    n_slots_ = n_slots;
  }

'''
    if marker not in s:
        raise SystemExit("Could not find Clear method in accumulator")
    s = s.replace(marker, marker + method, 1)

p.write_text(s)

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

start, end = find_function(s, "AddFlatIntermediateDirectionalValue")

new_body = """
inline bool AddFlatIntermediateDirectionalValue(const size_t direction,
                                                const VertexId i,
                                                const VertexId j,
                                                const Real value) {
  if (!g_ADGraph->useFlatIntermediateDirectionalBackend) {
    ++g_flat_intermediate_write_miss_count;
    RecordFlatIntermediateMissSample(
        g_flat_intermediate_write_miss_samples, i, j);
    return false;
  }

  const std::size_t slot =
      g_ADGraph->intermediateEdgeSlotRegistry.GetOrCreate(i, j);

  g_ADGraph->flatIntermediateDirectionalValues.EnsureSlotsPreserve(slot + 1);
  g_ADGraph->flatIntermediateDirectionalValues.Add(direction, slot, value);

  ++g_flat_intermediate_write_hit_count;
  return true;
}
"""

s = s[:start] + new_body + s[end:]
p.write_text(s)
PYEOF

python3 /tmp/repair_dynamic_slots.py

cat <<'EOF'

Installed dynamic intermediate flat slot allocation.

Run:
  ./run_had_quadra_nonzero_batch_directional_test.sh
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

Watch:
  flat intermediate read_hit/read_miss/write_hit/write_miss
  grad diff / obj diff

EOF
