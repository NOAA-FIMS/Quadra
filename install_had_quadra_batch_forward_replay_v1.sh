#!/usr/bin/env bash
set -euo pipefail

mkdir -p .quadra_patch_backups

target="core/had_quadra.hpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/had_quadra.hpp.batch_forward_replay.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/quadra_batch_forward_replay_patch.py <<'PYEOF'
from pathlib import Path

p = Path("core/had_quadra.hpp")
s = p.read_text()

if "PropagateDirectionalBatchForwardReplay" not in s:
    anchor = s.find("inline void PropagateAdjointDirectionalBatch()")
    if anchor < 0:
        raise SystemExit("Could not find PropagateAdjointDirectionalBatch")

    helper = """
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

    const Real lp =
        (left < g_ADGraph->vertices.size()) ? g_ADGraph->vertices[left].primal
                                            : Real(0.0);
    const Real rp =
        (right < g_ADGraph->vertices.size()) ? g_ADGraph->vertices[right].primal
                                             : Real(0.0);

    for (int k = 0; k < nDirections; ++k) {
      const size_t kk = static_cast<size_t>(k);
      const Real ld =
          (left < g_ADGraph->vertices.size() &&
           kk < g_ADGraph->vertices[left].dotBatch.size())
              ? g_ADGraph->vertices[left].dotBatch[kk]
              : Real(0.0);
      const Real rd =
          (right < g_ADGraph->vertices.size() &&
           kk < g_ADGraph->vertices[right].dotBatch.size())
              ? g_ADGraph->vertices[right].dotBatch[kk]
              : Real(0.0);

      switch (v.op) {
      case OpCode::Add:
        v.dotBatch[kk] = ld + rd;
        break;

      case OpCode::Sub:
        v.dotBatch[kk] = ld - rd;
        break;

      case OpCode::Mul:
        v.dotBatch[kk] = ld * rp + lp * rd;
        if (v.e1.to != vid) {
          v.e1.dwBatch[kk] = rd;
        }
        if (v.e2.to != vid) {
          v.e2.dwBatch[kk] = ld;
        }
        v.soWDotBatch[kk] = Real(0.0);
        break;

      case OpCode::Div:
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

      case OpCode::Sin:
        v.dotBatch[kk] = std::cos(lp) * ld;
        if (v.e1.to != vid) {
          v.e1.dwBatch[kk] = -std::sin(lp) * ld;
        }
        v.soWDotBatch[kk] = -std::cos(lp) * ld;
        break;

      case OpCode::Cos:
        v.dotBatch[kk] = -std::sin(lp) * ld;
        if (v.e1.to != vid) {
          v.e1.dwBatch[kk] = -std::cos(lp) * ld;
        }
        v.soWDotBatch[kk] = std::sin(lp) * ld;
        break;

      case OpCode::Neg:
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
    s = s[:anchor] + helper + s[anchor:]

idx = s.find("inline void PropagateAdjointDirectionalBatch()")
if idx < 0:
    raise SystemExit("Could not find batch reverse function")

guard = """  if (nDirections <= 0)
    return;"""
pos = s.find(guard, idx)
if pos < 0:
    raise SystemExit("Could not find nDirections guard")

after = pos + len(guard)
if "PropagateDirectionalBatchForwardReplay();" not in s[after:after+300]:
    s = s[:after] + "\n\n  PropagateDirectionalBatchForwardReplay();" + s[after:]

p.write_text(s)
PYEOF

python3 /tmp/quadra_batch_forward_replay_patch.py

cat <<'EOF'

Installed HAD Quadra batch forward replay v1.

Patched:
  core/had_quadra.hpp

Run:
  ./run_had_quadra_nonzero_batch_directional_test.sh
  ./run_had_quadra_directional_batch_propagation_test.sh
  ./run_exact_gradient_workspace_test.sh

EOF
