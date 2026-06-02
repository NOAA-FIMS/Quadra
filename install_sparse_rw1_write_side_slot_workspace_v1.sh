#!/usr/bin/env bash
set -euo pipefail

# install_sparse_rw1_write_side_slot_workspace_v1.sh
#
# Adds an experimental write-side directional slot workspace fast path for
# RW1/tridiagonal Hdot benchmarks.
#
# This is intentionally opt-in and benchmark-scoped:
#
#   - core/had_quadra.hpp gets generic optional batch slot hooks
#   - benchmark enables a tridiagonal slot map for random-effect vertices
#   - PushEdgeDotBatch writes to slots when a slot map is present
#   - falls back to BTree when no slot exists
#
# This tests the core idea:
#
#   soEdgesDotBatch[k][outer].Insert(inner, value)
#
# becomes:
#
#   slot_values[k][slot] += value
#
# for known Hdot pattern entries.

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

cp "$header" ".quadra_patch_backups/had_quadra.hpp.write_side_slot_workspace.$(date +%Y%m%d_%H%M%S).bak"
cp "$bench" ".quadra_patch_backups/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp.write_side_slot_workspace.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/quadra_write_side_slot_workspace_patch.py <<'PYEOF'
from pathlib import Path

h = Path("core/had_quadra.hpp")
s = h.read_text()

# Add storage to ADGraph.
if "batchSlotOuterInnerToSlot" not in s:
    old = """        // Indexed as soEdgesDotBatch[k][vertex] and selfSoEdgesDotBatch[k][vertex].
        std::vector<std::vector<BTree>> soEdgesDotBatch;
        std::vector<std::vector<Real>> selfSoEdgesDotBatch;"""
    new = """        // Indexed as soEdgesDotBatch[k][vertex] and selfSoEdgesDotBatch[k][vertex].
        std::vector<std::vector<BTree>> soEdgesDotBatch;
        std::vector<std::vector<Real>> selfSoEdgesDotBatch;

        // Optional write-side directional slot workspace.
        // batchSelfSlot[vertex] gives the direct slot for diagonal Hdot(vertex,vertex).
        // batchSlotOuterInnerToSlot[outer].Query(inner) gives off-diagonal slot,
        // or 0 if no slot is mapped. Slot ids are stored as slot+1 so 0 can mean absent.
        bool useBatchDirectionalSlotWorkspace = false;
        std::vector<int> batchSelfSlot;
        std::vector<BTree> batchSlotOuterInnerToSlot;
        std::vector<std::vector<Real>> batchDirectionalSlotValues;"""
    if old not in s:
        raise SystemExit("Could not find ADGraph batch storage block")
    s = s.replace(old, new, 1)

# Clear in graph clear if obvious location.
if "useBatchDirectionalSlotWorkspace = false;" not in s:
    old = """    soEdgesDot.clear();"""
    new = """    soEdgesDot.clear();
    useBatchDirectionalSlotWorkspace = false;
    batchSelfSlot.clear();
    batchSlotOuterInnerToSlot.clear();
    batchDirectionalSlotValues.clear();"""
    s = s.replace(old, new, 1)

# Add API helpers near batch helpers.
if "EnableBatchDirectionalSlotWorkspace" not in s:
    anchor = s.find("inline void ResizeDirectionalBatch")
    if anchor < 0:
        raise SystemExit("ResizeDirectionalBatch not found")
    helper = """
inline void DisableBatchDirectionalSlotWorkspace() {
  g_ADGraph->useBatchDirectionalSlotWorkspace = false;
  g_ADGraph->batchSelfSlot.clear();
  g_ADGraph->batchSlotOuterInnerToSlot.clear();
  g_ADGraph->batchDirectionalSlotValues.clear();
}

inline void EnableBatchDirectionalSlotWorkspace(const int nSlots) {
  const int nDirections = g_ADGraph->nBatchDirections;
  g_ADGraph->useBatchDirectionalSlotWorkspace = true;

  g_ADGraph->batchSelfSlot.assign(g_ADGraph->vertices.size(), -1);
  g_ADGraph->batchSlotOuterInnerToSlot.clear();
  g_ADGraph->batchSlotOuterInnerToSlot.resize(g_ADGraph->vertices.size());

  g_ADGraph->batchDirectionalSlotValues.assign(
      static_cast<size_t>(nDirections),
      std::vector<Real>(static_cast<size_t>(nSlots), Real(0.0)));
}

inline void ClearBatchDirectionalSlotValues() {
  if (!g_ADGraph->useBatchDirectionalSlotWorkspace) {
    return;
  }
  for (auto &values : g_ADGraph->batchDirectionalSlotValues) {
    std::fill(values.begin(), values.end(), Real(0.0));
  }
}

inline void SetBatchDirectionalSelfSlot(const VertexId vertex,
                                        const int slot) {
  if (g_ADGraph->batchSelfSlot.size() < g_ADGraph->vertices.size()) {
    g_ADGraph->batchSelfSlot.assign(g_ADGraph->vertices.size(), -1);
  }
  g_ADGraph->batchSelfSlot[vertex] = slot;
}

inline void SetBatchDirectionalOffdiagSlot(const VertexId i,
                                           const VertexId j,
                                           const int slot) {
  const VertexId outer = std::max(i, j);
  const VertexId inner = std::min(i, j);

  if (g_ADGraph->batchSlotOuterInnerToSlot.size() < g_ADGraph->vertices.size()) {
    g_ADGraph->batchSlotOuterInnerToSlot.resize(g_ADGraph->vertices.size());
  }

  // Store slot+1 so Query()==0 means absent.
  g_ADGraph->batchSlotOuterInnerToSlot[outer].Insert(inner, Real(slot + 1));
}

inline bool AddBatchDirectionalSlotValue(const size_t direction,
                                         const VertexId i,
                                         const VertexId j,
                                         const Real value) {
  if (!g_ADGraph->useBatchDirectionalSlotWorkspace) {
    return false;
  }
  if (direction >= g_ADGraph->batchDirectionalSlotValues.size()) {
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
    return false;
  }

  auto &values = g_ADGraph->batchDirectionalSlotValues[direction];
  if (slot >= static_cast<int>(values.size())) {
    return false;
  }

  values[static_cast<size_t>(slot)] += value;
  return true;
}

inline Real GetBatchDirectionalSlotValue(const size_t direction,
                                         const int slot) {
  return g_ADGraph->batchDirectionalSlotValues[direction][static_cast<size_t>(slot)];
}

"""
    s = s[:anchor] + helper + s[anchor:]

# Clear slots in ResizeDirectionalBatch after batch values assign if enabled.
if "if (g_ADGraph->useBatchDirectionalSlotWorkspace)" not in s:
    marker = """        for (int k = 0; k < nDirections; ++k)
        {
            g_ADGraph->soEdgesDotBatch[static_cast<size_t>(k)].resize(
                g_ADGraph->vertices.size());

            g_ADGraph->selfSoEdgesDotBatch[static_cast<size_t>(k)].assign(
                g_ADGraph->vertices.size(),
                Real(0.0));
        }"""
    repl = marker + """

        if (g_ADGraph->useBatchDirectionalSlotWorkspace)
        {
            const size_t nSlots =
                g_ADGraph->batchDirectionalSlotValues.empty()
                    ? 0
                    : g_ADGraph->batchDirectionalSlotValues.front().size();
            g_ADGraph->batchDirectionalSlotValues.assign(
                static_cast<size_t>(nDirections),
                std::vector<Real>(nSlots, Real(0.0)));
        }"""
    s = s.replace(marker, repl, 1)

# Patch PushEdgeDotBatch self/offdiag writes to direct slot first.
old = """    selfDots[foEdge.to] += Real(2.0) * valDot;"""
new = """    if (!AddBatchDirectionalSlotValue(direction, foEdge.to, foEdge.to,
                                      Real(2.0) * valDot)) {
      selfDots[foEdge.to] += Real(2.0) * valDot;
    }"""
s = s.replace(old, new, 1)

old = """    trees[outer].Insert(inner, valDot);"""
new = """    if (!AddBatchDirectionalSlotValue(direction, outer, inner, valDot)) {
      trees[outer].Insert(inner, valDot);
    }"""
s = s.replace(old, new, 1)

# Clear slot values before reverse sweep.
marker = """  PropagateDirectionalBatchForwardReplay();"""
if "ClearBatchDirectionalSlotValues();" not in s:
    s = s.replace(marker, marker + "\n  ClearBatchDirectionalSlotValues();", 1)

h.write_text(s)

# Patch benchmark adapter to enable slots and read trace from slots.
b = Path("benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp")
t = b.read_text()

# In adapter build(), after ResizeDirectionalBatch, install RW1 slot map.
old = """        workspace.ResizeDirectionalBatch(static_cast<std::size_t>(K));
    }"""
new = """        workspace.ResizeDirectionalBatch(static_cast<std::size_t>(K));

        // Experimental RW1/tridiagonal write-side directional slot workspace.
        // slots: 0..m-1 diag, m..2m-2 subdiag(i,i-1).
        workspace.HadWorkspace().Activate();
        had::EnableBatchDirectionalSlotWorkspace(2 * m - 1);

        for (int i = 0; i < m; ++i) {
            had::SetBatchDirectionalSelfSlot(
                random_vars[static_cast<std::size_t>(i)].varId,
                i);
            if (i > 0) {
                had::SetBatchDirectionalOffdiagSlot(
                    random_vars[static_cast<std::size_t>(i)].varId,
                    random_vars[static_cast<std::size_t>(i - 1)].varId,
                    m + i - 1);
            }
        }
    }"""
if old in t and "EnableBatchDirectionalSlotWorkspace" not in t:
    t = t.replace(old, new, 1)

# Replace exact_gradient trace assembly with direct slot reads if possible.
old = """        const auto assembled =
            workspace.AssembleExactGradient(
                joint_objective,
                logdet,
                joint_grad,
                std::forward<SelectedInverseAccessor>(selected_inverse),
                pattern);

        ExactGradientResult out;
        out.objective = assembled.objective;
        out.gradient = assembled.gradient;
        return out;"""
new = """        Eigen::VectorXd traces = Eigen::VectorXd::Zero(K);

        for (int k = 0; k < K; ++k) {
            double trace = 0.0;
            for (int i = 0; i < m; ++i) {
                trace += selected_inverse(i, i) *
                         had::GetBatchDirectionalSlotValue(
                             static_cast<size_t>(k), i);

                if (i > 0) {
                    trace += 2.0 * selected_inverse(i, i - 1) *
                             had::GetBatchDirectionalSlotValue(
                                 static_cast<size_t>(k), m + i - 1);
                }
            }
            traces[k] = trace;
        }

        ExactGradientResult out;
        out.objective = joint_objective + 0.5 * logdet;
        out.gradient = joint_grad + 0.5 * traces;
        return out;"""
if old in t:
    t = t.replace(old, new, 1)
else:
    print("WARNING: exact_gradient assembled block not found")

b.write_text(t)
PYEOF

python3 /tmp/quadra_write_side_slot_workspace_patch.py

cat <<'EOF'

Installed experimental RW1 write-side directional slot workspace.

Run:
  ./run_had_quadra_nonzero_batch_directional_test.sh
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

Notes:
  This is opt-in through the benchmark adapter and falls back to BTree
  for unmapped entries.

EOF
