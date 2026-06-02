#!/usr/bin/env bash
set -euo pipefail

# install_sparse_rw1_reuse_only_profile_mode_v1.sh
#
# Adds an optional reuse-only profile mode to:
#   benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp
#
# Usage:
#   ./build/benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse 10 --reuse-only
#
# In reuse-only mode:
#   - the rebuild reference loop is skipped during timing
#   - correctness is still checked once before timing if possible
#   - output focuses on the reuse buckets
#
# Purpose:
#   Time Profiler samples are currently mixed between exact_gradient_rebuild()
#   and ExactGradientWorkspace::PropagateDirectionalBatch(). This mode isolates
#   the optimized reuse path so the top symbols are cleaner.

mkdir -p .quadra_patch_backups

target="benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp.reuse_only.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/quadra_reuse_only_patch.py <<'PYEOF'
from pathlib import Path
import re

p = Path("benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp")
s = p.read_text()

# Add reuse_only parameter to run_case signature.
s = s.replace(
    "Row run_case(int m, int K, int reps) {",
    "Row run_case(int m, int K, int reps, bool reuse_only = false) {",
    1,
)

# Replace rebuild timing block with optional block.
old = """    const auto t0 = Clock::now();
    for (int r = 0; r < reps; ++r) {
        last_rebuild = exact_gradient_rebuild(
            m, K, theta, uhat, Hinv, joint_grad, joint_obj, logdet, factor);
    }
    const auto t1 = Clock::now();"""

new = """    const auto t0 = Clock::now();
    if (!reuse_only) {
        for (int r = 0; r < reps; ++r) {
            last_rebuild = exact_gradient_rebuild(
                m, K, theta, uhat, Hinv, joint_grad, joint_obj, logdet, factor);
        }
    } else {
        // In reuse-only profiling mode, compute the reference once outside the
        // timed loop so correctness is still checked without polluting samples.
        last_rebuild = exact_gradient_rebuild(
            m, K, theta, uhat, Hinv, joint_grad, joint_obj, logdet, factor);
    }
    const auto t1 = Clock::now();"""

if old not in s:
    raise SystemExit("Could not find rebuild timing block")
s = s.replace(old, new, 1)

# Adjust rebuild_ms assignment.
old = """    out.rebuild_ms = ms_between(t0, t1) / static_cast<double>(reps);"""
new = """    out.rebuild_ms = reuse_only
                         ? 0.0
                         : ms_between(t0, t1) / static_cast<double>(reps);"""
s = s.replace(old, new, 1)

# Parse --reuse-only in main.
old = """    int reps = 10;
    if (argc > 1) {
        reps = std::stoi(argv[1]);
    }"""

new = """    int reps = 10;
    bool reuse_only = false;

    if (argc > 1) {
        reps = std::stoi(argv[1]);
    }

    for (int i = 2; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--reuse-only") {
            reuse_only = true;
        }
    }"""
if old not in s:
    raise SystemExit("Could not find main arg parsing block")
s = s.replace(old, new, 1)

# Include string if needed.
if "#include <string>" not in s:
    s = s.replace("#include <iostream>", "#include <iostream>\n#include <string>", 1)

# Add mode print.
old = """    std::cout << "reps per case = " << reps << "\\n";
    std::cout << "K = " << K << "\\n\\n";"""
new = """    std::cout << "reps per case = " << reps << "\\n";
    std::cout << "K = " << K << "\\n";
    if (reuse_only) {
        std::cout << "mode = reuse-only profiling\\n";
    }
    std::cout << "\\n";"""
s = s.replace(old, new, 1)

# Pass reuse_only to run_case.
s = s.replace(
    "const Row r = run_case(m, K, reps);",
    "const Row r = run_case(m, K, reps, reuse_only);",
)

p.write_text(s)
PYEOF

python3 /tmp/quadra_reuse_only_patch.py

cat <<'EOF'

Installed reuse-only profile mode.

Run normally:
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

Build once, then profile reuse-only:
  ./build/benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse 10 --reuse-only

Profile command:
  xcrun xctrace record \
    --template "Time Profiler" \
    --output quadra_reuse_only.trace \
    --launch ./build/benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse \
    -- 10 --reuse-only

Then export:
  xcrun xctrace export \
    --input quadra_reuse_only.trace \
    --xpath '/trace-toc/run[@number="1"]/data/table[@schema="time-profile"]' \
    > quadra_reuse_only_profile.xml

EOF
