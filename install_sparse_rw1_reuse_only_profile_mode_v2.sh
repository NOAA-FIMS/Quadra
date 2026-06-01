#!/usr/bin/env bash
set -euo pipefail

mkdir -p .quadra_patch_backups

target="benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp.reuse_only_v2.$(date +%Y%m%d_%H%M%S).bak"

cat > /tmp/quadra_reuse_only_v2.py <<'PYEOF'
from pathlib import Path

p = Path("benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse.cpp")
s = p.read_text()

if "#include <string>" not in s:
    s = s.replace("#include <iostream>", "#include <iostream>\n#include <string>", 1)

s = s.replace(
    "Row run_case(int m, int K, int reps) {",
    "Row run_case(int m, int K, int reps, bool reuse_only = false) {",
    1,
)

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
        // One untimed reference evaluation for correctness without polluting
        // the profile with rebuild-path samples.
        last_rebuild = exact_gradient_rebuild(
            m, K, theta, uhat, Hinv, joint_grad, joint_obj, logdet, factor);
    }
    const auto t1 = Clock::now();"""

if old not in s:
    raise SystemExit("Could not find rebuild timing loop")
s = s.replace(old, new, 1)

s = s.replace(
    "    out.rebuild_ms = ms_between(t0, t1) / static_cast<double>(reps);",
    """    out.rebuild_ms = reuse_only
                         ? 0.0
                         : ms_between(t0, t1) / static_cast<double>(reps);""",
    1,
)

s = s.replace(
    "    out.speedup = out.rebuild_ms / out.reuse_ms;",
    "    out.speedup = reuse_only ? 0.0 : out.rebuild_ms / out.reuse_ms;",
    1,
)

old = """int main(int argc, char** argv) {
    int reps = 10;
    if (argc > 1) reps = std::stoi(argv[1]);"""

new = """int main(int argc, char** argv) {
    int reps = 10;
    bool reuse_only = false;

    if (argc > 1) reps = std::stoi(argv[1]);

    for (int i = 2; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--reuse-only") {
            reuse_only = true;
        }
    }"""

if old not in s:
    raise SystemExit("Could not find main arg parsing block")
s = s.replace(old, new, 1)

s = s.replace(
    """    std::cout << "reps per case = " << reps << "\\n";
    std::cout << "K = " << K << "\\n\\n";""",
    """    std::cout << "reps per case = " << reps << "\\n";
    std::cout << "K = " << K << "\\n";
    if (reuse_only) {
        std::cout << "mode = reuse-only profiling\\n";
    }
    std::cout << "\\n";""",
    1,
)

s = s.replace(
    "        const Row r = run_case(m, K, reps);",
    "        const Row r = run_case(m, K, reps, reuse_only);",
    1,
)

p.write_text(s)
PYEOF

python3 /tmp/quadra_reuse_only_v2.py

cat <<'EOF'

Installed reuse-only profile mode v2.

Run:
  ./run_sparse_rw1_exact_gradient_graph_reuse_benchmark.sh 10

Then profile only reuse:
  xcrun xctrace record \
    --template "Time Profiler" \
    --output quadra_reuse_only.trace \
    --launch ./build/benchmarks/benchmark_sparse_rw1_exact_gradient_graph_reuse \
    -- 10 --reuse-only

EOF
