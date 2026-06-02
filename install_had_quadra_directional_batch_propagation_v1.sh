#!/usr/bin/env bash
set -euo pipefail

# install_had_quadra_directional_batch_propagation_v1.sh
#
# Adds side-by-side:
#   had::PropagateAdjointDirectionalBatch()
#
# This does not replace PropagateAdjointDirectional().
# Correctness-first patch: scalar and batch paths must agree first.

mkdir -p tests .quadra_patch_backups

target="core/had_quadra.hpp"
if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

if ! grep -q "ResizeDirectionalBatch" "$target"; then
  echo "ERROR: batched directional scaffold is not installed yet."
  exit 1
fi

cp "$target" ".quadra_patch_backups/had_quadra.hpp.directional_batch_v1.$(date +%Y%m%d_%H%M%S).bak"

python3 - <<'PY'
from pathlib import Path
import re

p = Path("core/had_quadra.hpp")
s = p.read_text()

if "PropagateAdjointDirectionalBatch" in s:
    print("PropagateAdjointDirectionalBatch already exists; skipping header patch.")
    raise SystemExit(0)

m = re.search(r"\n\s*inline\s+void\s+PropagateAdjointDirectional\s*\(\s*\)\s*\{", s)
if not m:
    raise SystemExit("Could not find PropagateAdjointDirectional definition.")

brace = s.find("{", m.start())
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
    raise SystemExit("Could not find end of PropagateAdjointDirectional.")

batch_fn = r'''

    //==================================================
    // Batched directional reverse sweep.
    //
    // Correctness-first implementation:
    //   - preserves existing scalar PropagateAdjointDirectional()
    //   - mirrors scalar propagation over all active batch directions
    //   - writes batched second-order edge sensitivities into
    //     soEdgesDotBatch/selfSoEdgesDotBatch
    //==================================================

    inline void PropagateAdjointDirectionalBatch()
    {
        const int nDirections = g_ADGraph->nBatchDirections;

        if (nDirections <= 0)
        {
            return;
        }

        for (auto &v : g_ADGraph->vertices)
        {
            if (v.dotBatch.size() != static_cast<size_t>(nDirections))
            {
                v.dotBatch.resize(static_cast<size_t>(nDirections), Real(0.0));
            }

            v.wDotBatch.assign(static_cast<size_t>(nDirections), Real(0.0));
            v.soWDotBatch.assign(static_cast<size_t>(nDirections), Real(0.0));
            v.e1.dwBatch.assign(static_cast<size_t>(nDirections), Real(0.0));
            v.e2.dwBatch.assign(static_cast<size_t>(nDirections), Real(0.0));
        }

        g_ADGraph->soEdgesDotBatch.resize(static_cast<size_t>(nDirections));
        g_ADGraph->selfSoEdgesDotBatch.resize(static_cast<size_t>(nDirections));

        for (int k = 0; k < nDirections; ++k)
        {
            g_ADGraph->soEdgesDotBatch[static_cast<size_t>(k)].resize(
                g_ADGraph->vertices.size());

            for (auto &tree : g_ADGraph->soEdgesDotBatch[static_cast<size_t>(k)])
            {
                tree.Clear();
            }

            g_ADGraph->selfSoEdgesDotBatch[static_cast<size_t>(k)].assign(
                g_ADGraph->vertices.size(),
                Real(0.0));
        }

        for (auto it = g_ADGraph->vertices.rbegin();
             it != g_ADGraph->vertices.rend();
             ++it)
        {
            ADVertex &v = *it;

            for (int k = 0; k < nDirections; ++k)
            {
                const size_t kk = static_cast<size_t>(k);

                const Real vDot =
                    (v.dotBatch.size() > kk) ? v.dotBatch[kk] : Real(0.0);

                if (v.e1.parent != kInvalid && v.e1.parent < g_ADGraph->vertices.size())
                {
                    ADVertex &p = g_ADGraph->vertices[v.e1.parent];

                    const Real parentDot =
                        (p.dotBatch.size() > kk) ? p.dotBatch[kk] : Real(0.0);

                    const Real dw =
                        v.e1.w2 * parentDot + v.e1.w2_other * vDot;

                    v.e1.dwBatch[kk] = dw;

                    p.wDotBatch[kk] +=
                        v.wDotBatch[kk] * v.e1.w + v.w * dw;
                }

                if (v.e2.parent != kInvalid && v.e2.parent < g_ADGraph->vertices.size())
                {
                    ADVertex &p = g_ADGraph->vertices[v.e2.parent];

                    const Real parentDot =
                        (p.dotBatch.size() > kk) ? p.dotBatch[kk] : Real(0.0);

                    const Real dw =
                        v.e2.w2 * parentDot + v.e2.w2_other * vDot;

                    v.e2.dwBatch[kk] = dw;

                    p.wDotBatch[kk] +=
                        v.wDotBatch[kk] * v.e2.w + v.w * dw;
                }

                if (v.e1.parent != kInvalid && v.e1.parent < g_ADGraph->vertices.size())
                {
                    const VertexId child = v.varId;
                    const VertexId parent = v.e1.parent;
                    const VertexId outer = std::max(child, parent);
                    const VertexId inner = std::min(child, parent);

                    g_ADGraph->soEdgesDotBatch[kk][outer].Insert(
                        inner,
                        v.wDotBatch[kk] * v.e1.w + v.w * v.e1.dwBatch[kk]);
                }

                if (v.e2.parent != kInvalid && v.e2.parent < g_ADGraph->vertices.size())
                {
                    const VertexId child = v.varId;
                    const VertexId parent = v.e2.parent;
                    const VertexId outer = std::max(child, parent);
                    const VertexId inner = std::min(child, parent);

                    g_ADGraph->soEdgesDotBatch[kk][outer].Insert(
                        inner,
                        v.wDotBatch[kk] * v.e2.w + v.w * v.e2.dwBatch[kk]);
                }

                if (v.selfSecondOrder != Real(0.0) ||
                    v.selfSecondOrderDot != Real(0.0))
                {
                    g_ADGraph->selfSoEdgesDotBatch[kk][v.varId] +=
                        v.wDotBatch[kk] * v.selfSecondOrder +
                        v.w * v.selfSecondOrderDot;
                }
            }
        }
    }
'''

s = s[:end] + batch_fn + s[end:]
p.write_text(s)
PY

cat > tests/test_had_quadra_directional_batch_propagation.cpp <<'EOF'
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../core/had_quadra.hpp"

DECLARE_ADGRAPH()

namespace {

using had::AReal;
using had::Real;

AReal objective(const AReal& x, const AReal& y) {
    return exp(x * y) + x * x * y + sin(y);
}

Real scalar_directional_hxy(const Real x0,
                            const Real y0,
                            const Real dx,
                            const Real dy) {
    had::ADGraph graph;
    had::g_ADGraph = &graph;

    AReal x(x0);
    AReal y(y0);
    AReal f = objective(x, y);

    f.w = 1.0;
    had::PropagateAdjoint();

    x.dot = dx;
    y.dot = dy;

    had::PropagateAdjointDirectional();

    return had::GetAdjointDot(x, y);
}

Real batch_directional_hxy(const Real x0,
                           const Real y0,
                           const std::vector<Real>& dx,
                           const std::vector<Real>& dy,
                           const int k) {
    had::ADGraph graph;
    had::g_ADGraph = &graph;

    AReal x(x0);
    AReal y(y0);
    AReal f = objective(x, y);

    f.w = 1.0;
    had::PropagateAdjoint();

    had::ResizeDirectionalBatch(static_cast<int>(dx.size()));

    for (int j = 0; j < static_cast<int>(dx.size()); ++j) {
        had::SetARealDotBatch(x, j, dx[static_cast<size_t>(j)]);
        had::SetARealDotBatch(y, j, dy[static_cast<size_t>(j)]);
    }

    had::PropagateAdjointDirectionalBatch();

    return had::GetAdjointDotBatch(x, y, k);
}

}  // namespace

int main() {
    const Real x0 = 0.4;
    const Real y0 = -0.7;

    const std::vector<Real> dx{1.0, 0.0, 0.25};
    const std::vector<Real> dy{0.0, 1.0, -0.50};

    for (int k = 0; k < 3; ++k) {
        const Real scalar =
            scalar_directional_hxy(x0, y0, dx[static_cast<size_t>(k)], dy[static_cast<size_t>(k)]);

        const Real batch =
            batch_directional_hxy(x0, y0, dx, dy, k);

        const Real diff = std::abs(scalar - batch);

        std::cout << "direction " << k
                  << " scalar=" << scalar
                  << " batch=" << batch
                  << " diff=" << diff << "\n";

        if (diff > 1.0e-10) {
            throw std::runtime_error("batched directional propagation differs from scalar path");
        }
    }

    std::cout << "had_quadra directional batch propagation tests passed\n";
    return 0;
}
EOF

cat > run_had_quadra_directional_batch_propagation_test.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2}"

EIGEN_INCLUDE=""
if [[ -d external/Eigen ]]; then
  EIGEN_INCLUDE="-Iexternal/Eigen"
elif [[ -d core/eigen ]]; then
  EIGEN_INCLUDE="-Icore/eigen"
fi

mkdir -p build/tests

set -x
"${CXX}" ${CXXFLAGS} ${EIGEN_INCLUDE} -I. \
  tests/test_had_quadra_directional_batch_propagation.cpp \
  -o build/tests/test_had_quadra_directional_batch_propagation

./build/tests/test_had_quadra_directional_batch_propagation
EOF

chmod +x run_had_quadra_directional_batch_propagation_test.sh

cat <<'EOF'

Installed had_quadra directional batch propagation v1.

Patched:
  core/had_quadra.hpp

Files added:
  tests/test_had_quadra_directional_batch_propagation.cpp
  run_had_quadra_directional_batch_propagation_test.sh

Run:
  ./run_had_quadra_directional_batch_propagation_test.sh

Expected:
  had_quadra directional batch propagation tests passed

EOF
