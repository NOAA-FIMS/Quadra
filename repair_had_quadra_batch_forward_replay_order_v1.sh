#!/usr/bin/env bash
set -euo pipefail

# repair_had_quadra_batch_forward_replay_order_v1.sh
#
# Diagnosis:
#   PropagateDirectionalBatchForwardReplay() now produces non-zero dotBatch
#   on intermediate/output vertices, but Hdot readout remains zero.
#
# Likely cause:
#   PropagateAdjointDirectionalBatch() calls forward replay too early, then
#   later clears:
#     vertex.wDotBatch
#     vertex.soWDotBatch
#     vertex.e1.dwBatch
#     vertex.e2.dwBatch
#
# Fix:
#   Move PropagateDirectionalBatchForwardReplay() to after that reset block
#   and before the reverse vid loop.

mkdir -p .quadra_patch_backups

target="core/had_quadra.hpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/had_quadra.hpp.batch_forward_replay_order.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/repair_batch_forward_replay_order.py <<'PYEOF'
from pathlib import Path

p = Path("core/had_quadra.hpp")
s = p.read_text()

idx = s.find("inline void PropagateAdjointDirectionalBatch()")
if idx < 0:
    raise SystemExit("PropagateAdjointDirectionalBatch not found")

# Find function bounds.
brace = s.find("{", idx)
depth = 0
end = None
for i in range(brace, len(s)):
    if s[i] == "{":
        depth += 1
    elif s[i] == "}":
        depth -= 1
        if depth == 0:
            end = i
            break
if end is None:
    raise SystemExit("Could not find function end")

body = s[idx:end]

# Remove existing replay call(s) inside function body only.
body = body.replace("\n\n  PropagateDirectionalBatchForwardReplay();", "")
body = body.replace("\n  PropagateDirectionalBatchForwardReplay();", "")

# Insert after the per-vertex batch reset block, right before reverse sweep.
marker = "  for (VertexId vid = n_vertices - 1; vid > 0; --vid) {"
pos = body.find(marker)
if pos < 0:
    raise SystemExit("Could not find reverse vid loop marker")

insert = "  PropagateDirectionalBatchForwardReplay();\n\n"
body = body[:pos] + insert + body[pos:]

s = s[:idx] + body + s[end:]
p.write_text(s)
PYEOF

python3 /tmp/repair_batch_forward_replay_order.py

cat <<'EOF'

Repaired batch forward replay call order.

Run:
  ./run_had_quadra_nonzero_batch_directional_test.sh
  ./run_had_quadra_directional_batch_propagation_test.sh

EOF
