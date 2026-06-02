#!/usr/bin/env bash
set -euo pipefail

# repair_remove_broken_batch_reverse_timing_v1.sh
#
# Removes the broken intrusive timing instrumentation inserted into
# core/had_quadra.hpp by install_had_quadra_batch_reverse_timing_breakdown_v1.sh.
#
# It restores compilability while keeping the benchmark-side reverse timing
# that already exists.
#
# The earlier patch inserted timing snippets into the wrong scopes, causing:
#   - undeclared BatchReverseNowMs
#   - undeclared batchReverseSparseStart
#   - undeclared batchReverseDiagStart
#   - undeclared batchReverseCreateStart
#   - undeclared batchReverseAdjointStart
#
# This repair strips those snippets conservatively.

mkdir -p .quadra_patch_backups

header="core/had_quadra.hpp"
bench="benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp"

if [[ ! -f "$header" ]]; then
  echo "ERROR: missing $header"
  exit 1
fi

cp "$header" ".quadra_patch_backups/had_quadra.hpp.remove_broken_timing.$(date +%Y%m%d_%H%M%S).bak"
if [[ -f "$bench" ]]; then
  cp "$bench" ".quadra_patch_backups/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp.remove_broken_timing.$(date +%Y%m%d_%H%M%S).bak"
fi

cat > /tmp/remove_broken_batch_reverse_timing.py <<'PYEOF'
from pathlib import Path
import re

p = Path("core/had_quadra.hpp")
s = p.read_text()

# Remove helper declarations if present.
helper_pat = re.compile(
    r'\ninline bool g_batch_reverse_timing_enabled = false;.*?'
    r'inline void PrintBatchReverseTimingCounters\(\) \{.*?\n\}\n',
    re.S,
)
s = helper_pat.sub("\n", s)

# Remove setup start block.
s = re.sub(
    r'\n\s*const double batchReverseSetupStart =\n'
    r'\s*g_batch_reverse_timing_enabled \? BatchReverseNowMs\(\) : 0\.0;',
    '',
    s,
)

# Remove setup end block.
s = re.sub(
    r'\n\s*if \(g_batch_reverse_timing_enabled\) \{\n'
    r'\s*g_batch_reverse_setup_ms \+=\n'
    r'\s*BatchReverseNowMs\(\) - batchReverseSetupStart;\n'
    r'\s*\}\n',
    '\n',
    s,
)

# Remove sparse start.
s = re.sub(
    r'\n\s*const double batchReverseSparseStart =\n'
    r'\s*g_batch_reverse_timing_enabled \? BatchReverseNowMs\(\) : 0\.0;\n',
    '\n',
    s,
)

# Remove sparse end.
s = re.sub(
    r'\n\s*if \(g_batch_reverse_timing_enabled\) \{\n'
    r'\s*g_batch_reverse_sparse_edges_ms \+=\n'
    r'\s*BatchReverseNowMs\(\) - batchReverseSparseStart;\n'
    r'\s*\}\n',
    '\n',
    s,
)

# Remove diag start.
s = re.sub(
    r'\n\s*const double batchReverseDiagStart =\n'
    r'\s*g_batch_reverse_timing_enabled \? BatchReverseNowMs\(\) : 0\.0;\n',
    '\n',
    s,
)

# Remove diag end.
s = re.sub(
    r'\n\s*if \(g_batch_reverse_timing_enabled\) \{\n'
    r'\s*g_batch_reverse_diag_ms \+= BatchReverseNowMs\(\) - batchReverseDiagStart;\n'
    r'\s*\}\n',
    '\n',
    s,
)

# Remove create start.
s = re.sub(
    r'\n\s*const double batchReverseCreateStart =\n'
    r'\s*g_batch_reverse_timing_enabled \? BatchReverseNowMs\(\) : 0\.0;\n',
    '\n',
    s,
)

# Remove create end.
s = re.sub(
    r'\n\s*if \(g_batch_reverse_timing_enabled\) \{\n'
    r'\s*g_batch_reverse_create_ms \+=\n'
    r'\s*BatchReverseNowMs\(\) - batchReverseCreateStart;\n'
    r'\s*\}\n',
    '\n',
    s,
)

# Remove adjoint start.
s = re.sub(
    r'\n\s*const double batchReverseAdjointStart =\n'
    r'\s*g_batch_reverse_timing_enabled \? BatchReverseNowMs\(\) : 0\.0;\n',
    '\n',
    s,
)

# Remove adjoint end.
s = re.sub(
    r'\n\s*if \(g_batch_reverse_timing_enabled\) \{\n'
    r'\s*g_batch_reverse_adjoint_ms \+=\n'
    r'\s*BatchReverseNowMs\(\) - batchReverseAdjointStart;\n'
    r'\s*\}\n',
    '\n',
    s,
)

p.write_text(s)

# Remove benchmark calls to timing helpers if present.
b = Path("benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp")
if b.exists():
    t = b.read_text()
    t = re.sub(
        r'\n\s*const bool emit_reverse_timing = \(m == 500 && r == reps - 1\);\n'
        r'\s*if \(emit_reverse_timing\) \{\n'
        r'\s*had::ResetBatchReverseTimingCounters\(\);\n'
        r'\s*had::EnableBatchReverseTiming\(true\);\n'
        r'\s*\}\n',
        '\n',
        t,
    )
    t = re.sub(
        r'\n\s*if \(emit_reverse_timing\) \{\n'
        r'\s*had::EnableBatchReverseTiming\(false\);\n'
        r'\s*had::PrintBatchReverseTimingCounters\(\);\n'
        r'\s*\}\n',
        '\n',
        t,
    )
    b.write_text(t)
PYEOF

python3 /tmp/remove_broken_batch_reverse_timing.py

cat <<'EOF'

Removed broken intrusive batch reverse timing instrumentation.

Run:
  ./run_had_quadra_nonzero_batch_directional_test.sh
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

If both compile, we can add a safer targeted breakdown next by instrumenting a copied benchmark-local reverse routine rather than patching the production header in place.

EOF
