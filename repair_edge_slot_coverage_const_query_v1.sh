#!/usr/bin/env bash
set -euo pipefail

# repair_edge_slot_coverage_const_query_v1.sh
#
# Fixes:
#   error: 'this' argument to member function 'Query' has type 'const BTree',
#   but function is not marked const
#
# Cause:
#   BTree::Query() is non-const, so CountMappedKeys() cannot take the slot map
#   as const BTree&.
#
# Repair:
#   CountMappedKeys(BTree &slot_map_tree) const
#   and remove const from slot_map local in diagnostic.

mkdir -p .quadra_patch_backups

target="core/had_quadra.hpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/had_quadra.hpp.edge_slot_const_query_repair.$(date +%Y%m%d_%H%M%S).bak"

python3 - <<'PYEOF'
from pathlib import Path

p = Path("core/had_quadra.hpp")
s = p.read_text()

s = s.replace(
    "inline std::size_t CountMappedKeys(const BTree &slot_map_tree) const",
    "inline std::size_t CountMappedKeys(BTree &slot_map_tree) const",
)

s = s.replace(
    "    const auto &slot_map = g_ADGraph->batchSlotOuterInnerToSlot[outer];",
    "    auto &slot_map = g_ADGraph->batchSlotOuterInnerToSlot[outer];",
)

p.write_text(s)
PYEOF

cat <<'EOF'

Repaired edge-slot coverage diagnostic const Query issue.

Run:
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

EOF
