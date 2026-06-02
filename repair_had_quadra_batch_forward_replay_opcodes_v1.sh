#!/usr/bin/env bash
set -euo pipefail

# repair_had_quadra_batch_forward_replay_opcodes_v1.sh
#
# Repairs PropagateDirectionalBatchForwardReplay() to use the actual OpCode
# names in this had_quadra.hpp:
#
#   Subtract, Multiply, Divide, Negate
#   AddConstant, SubtractConstant, ConstantSubtract
#   MultiplyConstant, DivideConstant, ConstantDivide
#   Sqrt
#
# Also removes nonexistent Sin/Cos cases.

mkdir -p .quadra_patch_backups

target="core/had_quadra.hpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/had_quadra.hpp.repair_batch_forward_replay_opcodes.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/repair_batch_forward_replay_opcodes.py <<'PYEOF'
from pathlib import Path

p = Path("core/had_quadra.hpp")
s = p.read_text()

start = s.find("inline void PropagateDirectionalBatchForwardReplay()")
if start < 0:
    raise SystemExit("PropagateDirectionalBatchForwardReplay not found")

# Find function end.
brace = s.find("{", start)
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
    raise SystemExit("Could not find end of PropagateDirectionalBatchForwardReplay")

replacement = """
inline void PropagateDirectionalBatchForwardReplay() {
  const int nDirections = g_ADGraph->nBatchDirections;
  if (nDirections <= 0) {
    return;
  }

  const size_t batchSize = static_cast<size_t>(nDirections);
  const VertexId n_vertices =
      static_cast<VertexId>(g_ADGraph->vertices.size());

  for (VertexId vid = 0; vid < n_vertices; ++vid) {
    ADVertex &v = g_ADGraph->vertices[vid];

    if (v.dotBatch.size() != batchSize) {
      v.dotBatch.assign(batchSize, Real(0.0));
    }
    if (v.wDotBatch.size() != batchSize) {
      v.wDotBatch.assign(batchSize, Real(0.0));
    }
    if (v.soWDotBatch.size() != batchSize) {
      v.soWDotBatch.assign(batchSize, Real(0.0));
    }
    if (v.e1.dwBatch.size() != batchSize) {
      v.e1.dwBatch.assign(batchSize, Real(0.0));
    }
    if (v.e2.dwBatch.size() != batchSize) {
      v.e2.dwBatch.assign(batchSize, Real(0.0));
    }

    std::fill(v.wDotBatch.begin(), v.wDotBatch.end(), Real(0.0));
    std::fill(v.soWDotBatch.begin(), v.soWDotBatch.end(), Real(0.0));
    std::fill(v.e1.dwBatch.begin(), v.e1.dwBatch.end(), Real(0.0));
    std::fill(v.e2.dwBatch.begin(), v.e2.dwBatch.end(), Real(0.0));

    if (v.op == OpCode::Independent) {
      continue;
    }

    std::fill(v.dotBatch.begin(), v.dotBatch.end(), Real(0.0));

    const VertexId left = v.left;
    const VertexId right = v.right;

    const bool has_left = left < g_ADGraph->vertices.size();
    const bool has_right = right < g_ADGraph->vertices.size();

    const Real lp = has_left ? g_ADGraph->vertices[left].primal : Real(0.0);
    const Real rp = has_right ? g_ADGraph->vertices[right].primal : Real(0.0);
    const Real c = v.constant;

    for (int k = 0; k < nDirections; ++k) {
      const size_t kk = static_cast<size_t>(k);

      const Real ld =
          (has_left && kk < g_ADGraph->vertices[left].dotBatch.size())
              ? g_ADGraph->vertices[left].dotBatch[kk]
              : Real(0.0);

      const Real rd =
          (has_right && kk < g_ADGraph->vertices[right].dotBatch.size())
              ? g_ADGraph->vertices[right].dotBatch[kk]
              : Real(0.0);

      switch (v.op) {
      case OpCode::Add:
        v.dotBatch[kk] = ld + rd;
        break;

      case OpCode::AddConstant:
        v.dotBatch[kk] = ld;
        break;

      case OpCode::Subtract:
        v.dotBatch[kk] = ld - rd;
        break;

      case OpCode::SubtractConstant:
        v.dotBatch[kk] = ld;
        break;

      case OpCode::ConstantSubtract:
        v.dotBatch[kk] = -ld;
        break;

      case OpCode::Multiply:
        v.dotBatch[kk] = ld * rp + lp * rd;
        if (v.e1.to != vid) {
          v.e1.dwBatch[kk] = rd;
        }
        if (v.e2.to != vid) {
          v.e2.dwBatch[kk] = ld;
        }
        break;

      case OpCode::MultiplyConstant:
        v.dotBatch[kk] = c * ld;
        if (v.e1.to != vid) {
          v.e1.dwBatch[kk] = Real(0.0);
        }
        break;

      case OpCode::Divide:
        if (rp != Real(0.0)) {
          const Real inv = Real(1.0) / rp;
          const Real inv2 = inv * inv;
          const Real inv3 = inv2 * inv;
          v.dotBatch[kk] = (ld * rp - lp * rd) * inv2;

          if (v.e1.to != vid) {
            v.e1.dwBatch[kk] = -rd * inv2;
          }
          if (v.e2.to != vid) {
            v.e2.dwBatch[kk] =
                (-ld * inv2) + (Real(2.0) * lp * rd * inv3);
          }
          v.soWDotBatch[kk] = -rd * inv2;
        }
        break;

      case OpCode::DivideConstant:
        if (c != Real(0.0)) {
          v.dotBatch[kk] = ld / c;
        }
        break;

      case OpCode::ConstantDivide:
        if (lp != Real(0.0)) {
          const Real inv = Real(1.0) / lp;
          const Real inv2 = inv * inv;
          const Real inv3 = inv2 * inv;
          v.dotBatch[kk] = -c * ld * inv2;
          if (v.e1.to != vid) {
            v.e1.dwBatch[kk] = Real(2.0) * c * ld * inv3;
          }
          v.soWDotBatch[kk] = -Real(6.0) * c * ld * inv3 * inv;
        }
        break;

      case OpCode::Exp: {
        const Real ev = std::exp(lp);
        v.dotBatch[kk] = ev * ld;
        if (v.e1.to != vid) {
          v.e1.dwBatch[kk] = ev * ld;
        }
        v.soWDotBatch[kk] = ev * ld;
        break;
      }

      case OpCode::Log:
        if (lp != Real(0.0)) {
          const Real inv = Real(1.0) / lp;
          const Real inv2 = inv * inv;
          v.dotBatch[kk] = ld * inv;
          if (v.e1.to != vid) {
            v.e1.dwBatch[kk] = -ld * inv2;
          }
          v.soWDotBatch[kk] = Real(2.0) * ld * inv2 * inv;
        }
        break;

      case OpCode::Sqrt:
        if (lp > Real(0.0)) {
          const Real root = std::sqrt(lp);
          const Real inv_root = Real(1.0) / root;
          v.dotBatch[kk] = Real(0.5) * inv_root * ld;
          if (v.e1.to != vid) {
            v.e1.dwBatch[kk] =
                -Real(0.25) * inv_root / lp * ld;
          }
          v.soWDotBatch[kk] =
              Real(0.375) * inv_root / (lp * lp) * ld;
        }
        break;

      case OpCode::Negate:
        v.dotBatch[kk] = -ld;
        break;

      case OpCode::Independent:
      default:
        break;
      }
    }
  }
}
"""

s = s[:start] + replacement + s[end:]
p.write_text(s)
PYEOF

python3 /tmp/repair_batch_forward_replay_opcodes.py

cat <<'EOF'

Repaired batch forward replay OpCode names.

Run:
  ./run_had_quadra_nonzero_batch_directional_test.sh

EOF
