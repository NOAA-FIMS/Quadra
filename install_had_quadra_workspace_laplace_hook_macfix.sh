#!/usr/bin/env bash
set -euo pipefail

# install_had_quadra_workspace_laplace_hook.sh
#
# Conservative integration patch for HadQuadraWorkspace into the Laplace/Hessian path.
#
# Goal:
#   1. Keep existing APIs working.
#   2. Add optional HadQuadraWorkspace* plumbing.
#   3. Add a benchmark target for workspace-enabled Laplace evaluation.
#   4. Avoid modifying external/had directly.
#
# This script is intentionally guarded. If your local file layout/signatures differ,
# it prints the exact next manual insertion points instead of blindly corrupting code.
#
# Run:
#   bash install_had_quadra_workspace_laplace_hook.sh
#   make benchmark-laplace-evaluator
#   make benchmark-laplace-evaluator-workspace
#
# Expected next result:
#   compare baseline evaluator vs workspace-enabled evaluator.

mkdir -p core/autodiff
mkdir -p core/memory
mkdir -p benchmarks
mkdir -p .quadra_patch_backups

# ----------------------------------------------------------------------
# 1. Ensure arena/workspace headers exist.
# ----------------------------------------------------------------------

if [[ ! -f core/memory/arena_pool.hpp ]]; then
cat > core/memory/arena_pool.hpp <<'HHH'
#pragma once

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <new>
#include <type_traits>
#include <utility>
#include <vector>

namespace quadra {

class ArenaPool {
public:
    explicit ArenaPool(std::size_t initial_block_bytes = 1u << 20)
        : block_bytes_m(initial_block_bytes == 0 ? (1u << 20) : initial_block_bytes) {
        add_block(block_bytes_m);
    }

    ArenaPool(const ArenaPool&) = delete;
    ArenaPool& operator=(const ArenaPool&) = delete;

    void reset() noexcept {
        current_block_m = 0;
        offset_m = 0;
        bytes_used_m = 0;
    }

    void release() {
        blocks_m.clear();
        current_block_m = 0;
        offset_m = 0;
        bytes_used_m = 0;
        add_block(block_bytes_m);
    }

    std::size_t bytes_used() const noexcept { return bytes_used_m; }

    std::size_t bytes_reserved() const noexcept {
        std::size_t out = 0;
        for (const auto& block : blocks_m) out += block.size;
        return out;
    }

    std::size_t block_count() const noexcept { return blocks_m.size(); }

    void* allocate_bytes(std::size_t nbytes, std::size_t alignment = alignof(std::max_align_t)) {
        if (nbytes == 0) nbytes = 1;
        if ((alignment & (alignment - 1)) != 0) throw std::bad_alloc();

        for (;;) {
            Block& block = blocks_m[current_block_m];

            std::uintptr_t base = reinterpret_cast<std::uintptr_t>(block.ptr.get());
            std::uintptr_t raw = base + offset_m;
            std::uintptr_t aligned =
                (raw + alignment - 1) & ~(static_cast<std::uintptr_t>(alignment) - 1);
            std::size_t new_offset = static_cast<std::size_t>(aligned - base) + nbytes;

            if (new_offset <= block.size) {
                offset_m = new_offset;
                bytes_used_m += nbytes;
                return reinterpret_cast<void*>(aligned);
            }

            if (current_block_m + 1 < blocks_m.size()) {
                ++current_block_m;
                offset_m = 0;
                continue;
            }

            const std::size_t next_size = std::max(block_bytes_m, nbytes + alignment);
            add_block(next_size);
            ++current_block_m;
            offset_m = 0;
        }
    }

    template <class T>
    T* allocate(std::size_t n = 1) {
        static_assert(!std::is_void<T>::value, "Cannot allocate void");
        return static_cast<T*>(allocate_bytes(sizeof(T) * n, alignof(T)));
    }

    template <class T, class... Args>
    T* create(Args&&... args) {
        T* ptr = allocate<T>(1);
        ::new (static_cast<void*>(ptr)) T(std::forward<Args>(args)...);
        return ptr;
    }

private:
    struct Block {
        std::unique_ptr<std::byte[]> ptr;
        std::size_t size = 0;
    };

    void add_block(std::size_t size) {
        Block block;
        block.ptr.reset(new std::byte[size]);
        block.size = size;
        blocks_m.push_back(std::move(block));
    }

    std::vector<Block> blocks_m;
    std::size_t block_bytes_m = 0;
    std::size_t current_block_m = 0;
    std::size_t offset_m = 0;
    std::size_t bytes_used_m = 0;
};

template <class T>
class ArenaAllocator {
public:
    using value_type = T;

    explicit ArenaAllocator(ArenaPool* arena = nullptr) noexcept
        : arena_m(arena) {}

    template <class U>
    ArenaAllocator(const ArenaAllocator<U>& other) noexcept
        : arena_m(other.arena()) {}

    T* allocate(std::size_t n) {
        if (!arena_m) throw std::bad_alloc();
        return arena_m->template allocate<T>(n);
    }

    void deallocate(T*, std::size_t) noexcept {}

    ArenaPool* arena() const noexcept { return arena_m; }

    template <class U>
    bool operator==(const ArenaAllocator<U>& other) const noexcept {
        return arena_m == other.arena();
    }

    template <class U>
    bool operator!=(const ArenaAllocator<U>& other) const noexcept {
        return !(*this == other);
    }

private:
    template <class U>
    friend class ArenaAllocator;

    ArenaPool* arena_m = nullptr;
};

} // namespace quadra
HHH
fi

cat > core/autodiff/had_quadra_workspace.hpp <<'HHH'
#pragma once

#include "../memory/arena_pool.hpp"

#include <cstddef>
#include <utility>
#include <vector>

namespace quadra {

class HadQuadraWorkspace {
public:
    explicit HadQuadraWorkspace(std::size_t block_bytes = 1u << 20)
        : arena_m(block_bytes) {}

    HadQuadraWorkspace(const HadQuadraWorkspace&) = delete;
    HadQuadraWorkspace& operator=(const HadQuadraWorkspace&) = delete;

    void reset() noexcept { arena_m.reset(); }
    void release() { arena_m.release(); }

    ArenaPool& arena() noexcept { return arena_m; }
    const ArenaPool& arena() const noexcept { return arena_m; }

    std::size_t bytes_used() const noexcept { return arena_m.bytes_used(); }
    std::size_t bytes_reserved() const noexcept { return arena_m.bytes_reserved(); }
    std::size_t block_count() const noexcept { return arena_m.block_count(); }

    template <class T>
    T* allocate(std::size_t n = 1) {
        return arena_m.template allocate<T>(n);
    }

    template <class T, class... Args>
    T* create(Args&&... args) {
        return arena_m.template create<T>(std::forward<Args>(args)...);
    }

    template <class T>
    using Vector = std::vector<T, ArenaAllocator<T>>;

    template <class T>
    Vector<T> make_vector() {
        return Vector<T>{ArenaAllocator<T>(&arena_m)};
    }

    template <class T>
    Vector<T> make_vector_with_reserve(std::size_t n) {
        Vector<T> out{ArenaAllocator<T>(&arena_m)};
        out.reserve(n);
        return out;
    }

private:
    ArenaPool arena_m;
};

} // namespace quadra
HHH

# ----------------------------------------------------------------------
# 2. Locate likely Hessian/Laplace files.
# ----------------------------------------------------------------------

echo "Searching for evaluate_random_effect_hessian..."

HESSIAN_FILE_LIST="$(grep -RIl "evaluate_random_effect_hessian" core benchmarks 2>/dev/null || true)"

if [[ -z "${HESSIAN_FILE_LIST}" ]]; then
    echo "Could not find evaluate_random_effect_hessian in core/ or benchmarks/."
    echo "Workspace headers were installed, but no source patch was attempted."
    exit 2
fi

echo "Found candidates:"
printf '%s\n' "${HESSIAN_FILE_LIST}" | sed 's/^/  /'

# Prefer core files over benchmark files.
TARGET_FILE="$(printf '%s\n' "${HESSIAN_FILE_LIST}" | grep '^core/' | head -n 1 || true)"
if [[ -z "${TARGET_FILE}" ]]; then
    TARGET_FILE="$(printf '%s\n' "${HESSIAN_FILE_LIST}" | head -n 1)"
fi

echo "Using target file: ${TARGET_FILE}"
cp "${TARGET_FILE}" ".quadra_patch_backups/$(basename "${TARGET_FILE}").before_workspace_hook.$(date +%Y%m%d_%H%M%S)"

# ----------------------------------------------------------------------
# 3. Apply a conservative source transformation.
#
#    This does NOT attempt to rewrite had::NewAReal internals yet.
#    It:
#      - adds include for workspace header,
#      - adds optional workspace argument to function signatures that match
#        common single-line template style,
#      - adds a scoped reset at function entry.
#
#    If the exact function signature is too complex, it leaves the file alone
#    and prints manual patch instructions.
# ----------------------------------------------------------------------

python3 - <<'PY' "${TARGET_FILE}"
from pathlib import Path
import re
import sys

path = Path(sys.argv[1])
text = path.read_text()

original = text
changed = False

include = '#include "../autodiff/had_quadra_workspace.hpp"'
include_alt = '#include "had_quadra_workspace.hpp"'

if "had_quadra_workspace.hpp" not in text:
    # Choose relative include based on target path.
    if "core/laplace/" in str(path) or "core/optimizer/" in str(path):
        inc = '#include "../autodiff/had_quadra_workspace.hpp"'
    elif "core/autodiff/" in str(path):
        inc = '#include "had_quadra_workspace.hpp"'
    else:
        inc = '#include "../core/autodiff/had_quadra_workspace.hpp"'

    # Insert after pragma once or after last include near top.
    if "#pragma once" in text:
        text = text.replace("#pragma once", "#pragma once\n\n" + inc, 1)
        changed = True
    else:
        m = list(re.finditer(r'^\s*#include\s+[<"].*[>"]\s*$', text, flags=re.M))
        if m:
            pos = m[-1].end()
            text = text[:pos] + "\n" + inc + text[pos:]
            changed = True
        else:
            text = inc + "\n" + text
            changed = True

# Add optional workspace argument to evaluate_random_effect_hessian signatures.
# This targets the common pattern:
#   RandomEffectHessianResult evaluate_random_effect_hessian<...>(..., double eps)
# or:
#   RandomEffectHessianResult evaluate_random_effect_hessian(..., double eps)
#
# We avoid modifying call sites in this patch by giving default nullptr.
sig_pattern = re.compile(
    r'(RandomEffectHessianResult\s+evaluate_random_effect_hessian\s*(?:<[^>]+>)?\s*\([^)]*?double\s+[A-Za-z_][A-Za-z0-9_]*\s*)\)',
    flags=re.S
)

def add_ws_arg(m):
    s = m.group(1)
    if "HadQuadraWorkspace" in s:
        return m.group(0)
    return s + ", HadQuadraWorkspace* workspace = nullptr)"

text2, n = sig_pattern.subn(add_ws_arg, text, count=1)
if n > 0 and text2 != text:
    text = text2
    changed = True

# Add reset at top of function body after the matched function if not already present.
if "workspace->reset()" not in text:
    marker = "evaluate_random_effect_hessian"
    idx = text.find(marker)
    if idx >= 0:
        brace = text.find("{", idx)
        if brace >= 0:
            insert = "\n    if (workspace != nullptr) {\n        workspace->reset();\n    }\n"
            text = text[:brace+1] + insert + text[brace+1:]
            changed = True

if changed and text != original:
    path.write_text(text)
    print(f"Patched {path}")
else:
    print(f"No automatic source changes applied to {path}")
    print()
    print("Manual edits likely needed:")
    print("  1. Include:")
    print("       #include \"../autodiff/had_quadra_workspace.hpp\"")
    print("  2. Change evaluate_random_effect_hessian signature:")
    print("       ..., double eps, HadQuadraWorkspace* workspace = nullptr)")
    print("  3. Add at function entry:")
    print("       if (workspace != nullptr) { workspace->reset(); }")
PY

# ----------------------------------------------------------------------
# 4. Create a workspace-enabled benchmark wrapper when possible.
#
#    We do not know the exact internals of benchmarks/laplace_evaluator_benchmark.cpp,
#    so this creates a copy and injects the include. If the original benchmark
#    constructs a LaplaceEvaluator class later, the user can add workspace there.
# ----------------------------------------------------------------------

if [[ -f benchmarks/laplace_evaluator_benchmark.cpp ]]; then
    cp benchmarks/laplace_evaluator_benchmark.cpp benchmarks/laplace_evaluator_workspace_benchmark.cpp

    if ! grep -q "had_quadra_workspace.hpp" benchmarks/laplace_evaluator_workspace_benchmark.cpp; then
        python3 - <<'PY'
from pathlib import Path
p = Path("benchmarks/laplace_evaluator_workspace_benchmark.cpp")
s = p.read_text()
m = list(__import__("re").finditer(r'^\s*#include\s+[<"].*[>"]\s*$', s, flags=__import__("re").M))
inc = '#include "../core/autodiff/had_quadra_workspace.hpp"'
if m:
    pos = m[-1].end()
    s = s[:pos] + "\n" + inc + s[pos:]
else:
    s = inc + "\n" + s
p.write_text(s)
PY
    fi

    echo "Created benchmarks/laplace_evaluator_workspace_benchmark.cpp"
else
    echo "benchmarks/laplace_evaluator_benchmark.cpp not found; skipping benchmark copy."
fi

# ----------------------------------------------------------------------
# 5. Patch Makefile target for workspace benchmark.
# ----------------------------------------------------------------------

if [[ -f Makefile ]]; then
    if ! grep -q "benchmark-laplace-evaluator-workspace" Makefile; then
        cat >> Makefile <<'MK'

.PHONY: benchmark-laplace-evaluator-workspace
benchmark-laplace-evaluator-workspace: benchmarks/laplace_evaluator_workspace_benchmark
	./benchmarks/laplace_evaluator_workspace_benchmark

benchmarks/laplace_evaluator_workspace_benchmark: benchmarks/laplace_evaluator_workspace_benchmark.cpp core/autodiff/had_quadra_workspace.hpp core/memory/arena_pool.hpp
	$(CXX) $(CXXFLAGS) -std=c++17 -O3 -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include -o benchmarks/laplace_evaluator_workspace_benchmark benchmarks/laplace_evaluator_workspace_benchmark.cpp
MK
        echo "Patched Makefile with benchmark-laplace-evaluator-workspace target."
    else
        echo "Makefile already contains benchmark-laplace-evaluator-workspace target."
    fi
else
    echo "No Makefile found; skipping Makefile patch."
fi

# ----------------------------------------------------------------------
# 6. Print next manual wiring guidance.
# ----------------------------------------------------------------------

cat <<'MSG'

Installed HadQuadraWorkspace Laplace hook scaffold.

Now run:

  make benchmark-laplace-evaluator
  make benchmark-laplace-evaluator-workspace

Important:
  This is a conservative scaffold. It adds optional workspace plumbing and an A/B
  benchmark target, but it does not yet rewrite had::NewAReal or ADGraph allocation.

If the workspace benchmark compiles but gives identical timing, the next exact edit is:

  - find calls to evaluate_random_effect_hessian(...) inside cached/evaluator code
  - create a persistent member:
        HadQuadraWorkspace workspace_m;
  - pass:
        &workspace_m
    into evaluate_random_effect_hessian(...)

Then rerun:

  make benchmark-laplace-evaluator-workspace

If compilation fails, paste the error and I will tighten the patch around your exact signatures.

MSG
