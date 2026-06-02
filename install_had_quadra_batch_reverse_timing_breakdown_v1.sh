#!/usr/bin/env bash
set -euo pipefail

mkdir -p .quadra_patch_backups

header="core/had_quadra.hpp"
bench="benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp"

if [[ ! -f "$header" ]]; then echo "ERROR: missing $header"; exit 1; fi
if [[ ! -f "$bench" ]]; then echo "ERROR: missing $bench"; exit 1; fi

cp "$header" ".quadra_patch_backups/had_quadra.hpp.batch_reverse_timing.$(date +%Y%m%d_%H%M%S).bak"
cp "$bench" ".quadra_patch_backups/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp.batch_reverse_timing.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/batch_reverse_timing_patch.py <<'PYEOF'
from pathlib import Path

h = Path("core/had_quadra.hpp")
s = h.read_text()

if "#include <chrono>" not in s:
    if "#include <iostream>" in s:
        s = s.replace("#include <iostream>", "#include <iostream>\n#include <chrono>", 1)
    elif "#include <cmath>" in s:
        s = s.replace("#include <cmath>", "#include <cmath>\n#include <chrono>", 1)
    else:
        s = "#include <chrono>\n" + s

if "g_batch_reverse_timing_enabled" not in s:
    anchor = s.find("inline std::uint64_t g_batch_query_count")
    if anchor < 0:
        raise SystemExit("Could not find batch counter anchor")

    block = """
inline bool g_batch_reverse_timing_enabled = false;
inline double g_batch_reverse_setup_ms = 0.0;
inline double g_batch_reverse_sparse_edges_ms = 0.0;
inline double g_batch_reverse_diag_ms = 0.0;
inline double g_batch_reverse_create_ms = 0.0;
inline double g_batch_reverse_adjoint_ms = 0.0;

inline double BatchReverseNowMs() {
  using clock = std::chrono::steady_clock;
  return std::chrono::duration<double, std::milli>(
      clock::now().time_since_epoch()).count();
}

inline void EnableBatchReverseTiming(const bool enabled) {
  g_batch_reverse_timing_enabled = enabled;
}

inline void ResetBatchReverseTimingCounters() {
  g_batch_reverse_setup_ms = 0.0;
  g_batch_reverse_sparse_edges_ms = 0.0;
  g_batch_reverse_diag_ms = 0.0;
  g_batch_reverse_create_ms = 0.0;
  g_batch_reverse_adjoint_ms = 0.0;
}

inline void PrintBatchReverseTimingCounters() {
  std::cerr << "HAD batch reverse timing ms:"
            << " setup=" << g_batch_reverse_setup_ms
            << " sparse_edges=" << g_batch_reverse_sparse_edges_ms
            << " diag=" << g_batch_reverse_diag_ms
            << " create=" << g_batch_reverse_create_ms
            << " adjoint=" << g_batch_reverse_adjoint_ms
            << " total="
            << (g_batch_reverse_setup_ms +
                g_batch_reverse_sparse_edges_ms +
                g_batch_reverse_diag_ms +
                g_batch_reverse_create_ms +
                g_batch_reverse_adjoint_ms)
            << "\\n";
}

"""
    s = s[:anchor] + block + s[anchor:]

# Add a setup timer after guard if not present.
old = """  if (nDirections <= 0)
    return;"""
new = """  if (nDirections <= 0)
    return;

  const double batchReverseSetupStart =
      g_batch_reverse_timing_enabled ? BatchReverseNowMs() : 0.0;"""
if "batchReverseSetupStart" not in s:
    s = s.replace(old, new, 1)

# End setup before reverse loop.
idx = s.find("inline void PropagateAdjointDirectionalBatch()")
marker = "  for (VertexId vid = n_vertices - 1; vid > 0; --vid) {"
pos = s.find(marker, idx)
if pos < 0:
    raise SystemExit("Could not find reverse loop")
if "g_batch_reverse_setup_ms +=\n        BatchReverseNowMs() - batchReverseSetupStart" not in s[idx:pos]:
    insert = """  if (g_batch_reverse_timing_enabled) {
    g_batch_reverse_setup_ms +=
        BatchReverseNowMs() - batchReverseSetupStart;
  }

"""
    s = s[:pos] + insert + s[pos:]

# Sparse edge timing.
old = """    for (auto it = btree.nodes.begin(); it != btree.nodes.end(); ++it) {
      ADEdge soEdge(it->key, it->val);"""
new = """    const double batchReverseSparseStart =
        g_batch_reverse_timing_enabled ? BatchReverseNowMs() : 0.0;

    for (auto it = btree.nodes.begin(); it != btree.nodes.end(); ++it) {
      ADEdge soEdge(it->key, it->val);"""
if "batchReverseSparseStart" not in s:
    s = s.replace(old, new, 1)

old = """    }

    const Real S = g_ADGraph->selfSoEdges[vid];"""
new = """    }

    if (g_batch_reverse_timing_enabled) {
      g_batch_reverse_sparse_edges_ms +=
          BatchReverseNowMs() - batchReverseSparseStart;
    }

    const double batchReverseDiagStart =
        g_batch_reverse_timing_enabled ? BatchReverseNowMs() : 0.0;

    const Real S = g_ADGraph->selfSoEdges[vid];"""
if "batchReverseDiagStart" not in s:
    s = s.replace(old, new, 1)

old = """    const Real a = vertex.w;"""
new = """    if (g_batch_reverse_timing_enabled) {
      g_batch_reverse_diag_ms += BatchReverseNowMs() - batchReverseDiagStart;
    }

    const double batchReverseCreateStart =
        g_batch_reverse_timing_enabled ? BatchReverseNowMs() : 0.0;

    const Real a = vertex.w;"""
if "batchReverseCreateStart" not in s:
    s = s.replace(old, new, 1)

old = """    if (a != Real(0.0)) {
      vertex.w = Real(0.0);"""
new = """    if (g_batch_reverse_timing_enabled) {
      g_batch_reverse_create_ms +=
          BatchReverseNowMs() - batchReverseCreateStart;
    }

    const double batchReverseAdjointStart =
        g_batch_reverse_timing_enabled ? BatchReverseNowMs() : 0.0;

    if (a != Real(0.0)) {
      vertex.w = Real(0.0);"""
if "batchReverseAdjointStart" not in s:
    s = s.replace(old, new, 1)

# Insert adjoint end before the loop closes. Use targeted tail.
tail = """    }
  }
}"""
insert_tail = """    }

    if (g_batch_reverse_timing_enabled) {
      g_batch_reverse_adjoint_ms +=
          BatchReverseNowMs() - batchReverseAdjointStart;
    }
  }
}"""
idx = s.find("inline void PropagateAdjointDirectionalBatch()")
end_idx = s.find("\n}\n", idx)
if "g_batch_reverse_adjoint_ms +=" not in s[idx:idx+12000]:
    last = s.rfind(tail, idx, idx+20000)
    if last >= 0:
        s = s[:last] + insert_tail + s[last+len(tail):]
    else:
        print("WARNING: could not insert adjoint bucket end")

h.write_text(s)

b = Path("benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp")
t = b.read_text()

old = """        had::ResetBatchDirectionalCounters();
        workspace.propagate_directional_batch();
        const auto rb3 = Clock::now();"""
new = """        had::ResetBatchDirectionalCounters();
        const bool emit_reverse_timing = (m == 500 && r == reps - 1);
        if (emit_reverse_timing) {
            had::ResetBatchReverseTimingCounters();
            had::EnableBatchReverseTiming(true);
        }

        workspace.propagate_directional_batch();

        if (emit_reverse_timing) {
            had::EnableBatchReverseTiming(false);
            had::PrintBatchReverseTimingCounters();
        }

        const auto rb3 = Clock::now();"""
if old not in t:
    raise SystemExit("Could not find benchmark reverse call block")
t = t.replace(old, new, 1)

b.write_text(t)
PYEOF

python3 /tmp/batch_reverse_timing_patch.py

cat <<'EOF'

Installed HAD batch reverse timing breakdown.

Run:
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

Look for stderr line:
  HAD batch reverse timing ms: setup=... sparse_edges=... diag=... create=... adjoint=...

EOF
