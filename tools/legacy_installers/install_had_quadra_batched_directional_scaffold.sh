#!/usr/bin/env bash
set -euo pipefail

mkdir -p tests .quadra_patch_backups
target="core/had_quadra.hpp"

if [[ ! -f "$target" ]]; then
  echo "ERROR: missing $target"
  exit 1
fi

cp "$target" ".quadra_patch_backups/had_quadra.hpp.batched_directional_scaffold.$(date +%Y%m%d_%H%M%S).bak"

python3 - <<'PY'
from pathlib import Path
import re

path = Path("core/had_quadra.hpp")
text = path.read_text()

if "ResizeDirectionalBatch" in text:
    print("had_quadra.hpp already contains batched directional scaffold; skipping.")
    raise SystemExit(0)

if "#include <vector>" not in text:
    matches = list(re.finditer(r"^#include .*$", text, flags=re.MULTILINE))
    if matches:
        pos = matches[-1].end()
        text = text[:pos] + "\n#include <vector>" + text[pos:]
    else:
        text = "#include <vector>\n" + text

def replace_once(haystack, needle, repl, label):
    if needle not in haystack:
        raise SystemExit(f"Could not find patch location: {label}")
    return haystack.replace(needle, repl, 1)

text = replace_once(
    text,
    "Real dot = Real(0.0);",
    '''Real dot = Real(0.0);

        // Batched first-order directional tangents for future multi-direction
        // Hdot propagation. The scalar dot remains the active production path.
        std::vector<Real> dotBatch;''',
    "AReal dot",
)

text = replace_once(
    text,
    "Real dw; // directional derivative of edge weight",
    '''Real dw; // directional derivative of edge weight

        // Batched directional derivative of edge weight.
        std::vector<Real> dwBatch;''',
    "ADEdge dw",
)

text = replace_once(
    text,
    '''        // Optional directional tangent associated with this vertex.
        Real dot;''',
    '''        // Optional directional tangent associated with this vertex.
        Real dot;

        // Batched directional tangent associated with this vertex.
        std::vector<Real> dotBatch;

        // Batched directional derivative of first-order adjoint weight.
        std::vector<Real> wDotBatch;

        // Batched directional derivative of soW along seeded primal tangents.
        std::vector<Real> soWDotBatch;''',
    "ADVertex dot",
)

text = replace_once(
    text,
    '''        std::vector<BTree> soEdgesDot;
        std::vector<Real> selfSoEdgesDot;''',
    '''        std::vector<BTree> soEdgesDot;
        std::vector<Real> selfSoEdgesDot;

        // Number of active batched directional tangents.
        int nBatchDirections = 0;

        // Batched directional Hessian edge storage.
        // Indexed as soEdgesDotBatch[k][vertex] and selfSoEdgesDotBatch[k][vertex].
        std::vector<std::vector<BTree>> soEdgesDotBatch;
        std::vector<std::vector<Real>> selfSoEdgesDotBatch;''',
    "ADGraph dot edge storage",
)

anchor = "    inline void PropagateAdjointDirectional()"
if anchor not in text:
    raise SystemExit("Could not find PropagateAdjointDirectional anchor.")

helpers = r'''
    //==================================================
    // Batched directional propagation scaffold.
    //
    // These helpers add graph storage and public accessors for multiple
    // simultaneous directional tangents. They intentionally do not replace
    // PropagateAdjointDirectional() yet.
    //==================================================

    inline void ResizeDirectionalBatch(const int nDirections)
    {
        if (nDirections < 0)
        {
            throw std::invalid_argument(
                "ResizeDirectionalBatch: nDirections must be nonnegative");
        }

        g_ADGraph->nBatchDirections = nDirections;

        for (auto &v : g_ADGraph->vertices)
        {
            v.dotBatch.assign(static_cast<size_t>(nDirections), Real(0.0));
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

            g_ADGraph->selfSoEdgesDotBatch[static_cast<size_t>(k)].assign(
                g_ADGraph->vertices.size(),
                Real(0.0));
        }
    }

    inline void ClearDirectionalBatch()
    {
        const int nDirections = g_ADGraph->nBatchDirections;

        for (auto &v : g_ADGraph->vertices)
        {
            v.dotBatch.assign(static_cast<size_t>(nDirections), Real(0.0));
            v.wDotBatch.assign(static_cast<size_t>(nDirections), Real(0.0));
            v.soWDotBatch.assign(static_cast<size_t>(nDirections), Real(0.0));
            v.e1.dwBatch.assign(static_cast<size_t>(nDirections), Real(0.0));
            v.e2.dwBatch.assign(static_cast<size_t>(nDirections), Real(0.0));
        }

        for (int k = 0; k < nDirections; ++k)
        {
            auto &trees = g_ADGraph->soEdgesDotBatch[static_cast<size_t>(k)];
            for (auto &tree : trees)
            {
                tree.Clear();
            }

            g_ADGraph->selfSoEdgesDotBatch[static_cast<size_t>(k)].assign(
                g_ADGraph->vertices.size(),
                Real(0.0));
        }
    }

    inline void CheckDirectionalBatchIndex(const int k)
    {
        if (k < 0 || k >= g_ADGraph->nBatchDirections)
        {
            throw std::out_of_range("directional batch index out of range");
        }
    }

    inline void SetVertexDotBatch(const VertexId vertexId,
                                  const int k,
                                  const Real value)
    {
        CheckDirectionalBatchIndex(k);

        if (vertexId >= g_ADGraph->vertices.size())
        {
            throw std::out_of_range("SetVertexDotBatch: vertexId out of range");
        }

        auto &v = g_ADGraph->vertices[vertexId];
        if (v.dotBatch.size() !=
            static_cast<size_t>(g_ADGraph->nBatchDirections))
        {
            v.dotBatch.assign(
                static_cast<size_t>(g_ADGraph->nBatchDirections),
                Real(0.0));
        }

        v.dotBatch[static_cast<size_t>(k)] = value;
    }

    inline Real GetVertexDotBatch(const VertexId vertexId,
                                  const int k)
    {
        CheckDirectionalBatchIndex(k);

        if (vertexId >= g_ADGraph->vertices.size())
        {
            throw std::out_of_range("GetVertexDotBatch: vertexId out of range");
        }

        const auto &v = g_ADGraph->vertices[vertexId];

        if (v.dotBatch.size() <= static_cast<size_t>(k))
        {
            return Real(0.0);
        }

        return v.dotBatch[static_cast<size_t>(k)];
    }

    inline void SetARealDotBatch(AReal &x,
                                 const int k,
                                 const Real value)
    {
        SetVertexDotBatch(x.varId, k, value);

        if (x.dotBatch.size() !=
            static_cast<size_t>(g_ADGraph->nBatchDirections))
        {
            x.dotBatch.assign(
                static_cast<size_t>(g_ADGraph->nBatchDirections),
                Real(0.0));
        }

        x.dotBatch[static_cast<size_t>(k)] = value;
    }

    inline Real GetAdjointDotBatch(const AReal &i,
                                   const AReal &j,
                                   const int k)
    {
        CheckDirectionalBatchIndex(k);

        if (g_ADGraph->soEdgesDotBatch.size() <= static_cast<size_t>(k) ||
            g_ADGraph->selfSoEdgesDotBatch.size() <= static_cast<size_t>(k))
        {
            return Real(0.0);
        }

        if (i.varId == j.varId)
        {
            const auto &diag =
                g_ADGraph->selfSoEdgesDotBatch[static_cast<size_t>(k)];

            if (diag.size() <= i.varId)
            {
                return Real(0.0);
            }

            return diag[i.varId];
        }

        const auto &trees =
            g_ADGraph->soEdgesDotBatch[static_cast<size_t>(k)];

        const VertexId outer = std::max(i.varId, j.varId);
        const VertexId inner = std::min(i.varId, j.varId);

        if (trees.size() <= outer)
        {
            return Real(0.0);
        }

        return trees[outer].Query(inner);
    }

'''
text = text.replace(anchor, helpers + anchor, 1)
path.write_text(text)
PY

cat > tests/test_had_quadra_batched_directional_scaffold.cpp <<'EOF'
#include <iostream>
#include <stdexcept>

#include "../core/had_quadra.hpp"

DECLARE_ADGRAPH()

int main()
{
    had::ADGraph graph;
    had::ADScope scope(graph);

    had::AReal x(2.0);
    had::AReal y(3.0);

    had::ResizeDirectionalBatch(3);

    had::SetARealDotBatch(x, 0, 1.0);
    had::SetARealDotBatch(x, 1, 2.0);
    had::SetARealDotBatch(y, 2, 4.0);

    if (had::GetVertexDotBatch(x.varId, 0) != 1.0)
    {
        throw std::runtime_error("x direction 0 was not stored");
    }

    if (had::GetVertexDotBatch(x.varId, 1) != 2.0)
    {
        throw std::runtime_error("x direction 1 was not stored");
    }

    if (had::GetVertexDotBatch(y.varId, 2) != 4.0)
    {
        throw std::runtime_error("y direction 2 was not stored");
    }

    had::ClearDirectionalBatch();

    if (had::GetVertexDotBatch(x.varId, 0) != 0.0)
    {
        throw std::runtime_error("clear did not reset x direction 0");
    }

    if (had::GetAdjointDotBatch(x, y, 0) != 0.0)
    {
        throw std::runtime_error("empty batched adjoint dot should be zero");
    }

    bool threw = false;
    try
    {
        had::SetARealDotBatch(x, 4, 1.0);
    }
    catch (const std::out_of_range &)
    {
        threw = true;
    }

    if (!threw)
    {
        throw std::runtime_error("out-of-range batch index did not throw");
    }

    std::cout << "had_quadra batched directional scaffold tests passed\n";
    return 0;
}
EOF

cat > run_had_quadra_batched_directional_scaffold_test.sh <<'EOF'
#!/usr/bin/env bash
set -euo pipefail

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2}"

EIGEN_INCLUDE=""
if [[ -d external/Eigen ]]; then
  EIGEN_INCLUDE="-Iexternal/Eigen"
elif [[ -d core/eigen ]]; then
  EIGEN_INCLUDE="-Icore/eigen"
elif [[ -d eigen ]]; then
  EIGEN_INCLUDE="-Ieigen"
elif [[ -d Eigen ]]; then
  EIGEN_INCLUDE="-I."
fi

LBFGS_INCLUDE=""
if [[ -d external/LBFGSpp/include ]]; then
  LBFGS_INCLUDE="-Iexternal/LBFGSpp/include"
elif [[ -d external/LBFGSpp ]]; then
  LBFGS_INCLUDE="-Iexternal/LBFGSpp"
fi

mkdir -p build/tests

set -x
"${CXX}" ${CXXFLAGS} ${EIGEN_INCLUDE} ${LBFGS_INCLUDE} -I.   tests/test_had_quadra_batched_directional_scaffold.cpp   -o build/tests/test_had_quadra_batched_directional_scaffold

./build/tests/test_had_quadra_batched_directional_scaffold
EOF

chmod +x run_had_quadra_batched_directional_scaffold_test.sh

cat <<'EOF'

Installed had_quadra batched directional scaffold.

Patched:
  core/had_quadra.hpp

Files added:
  tests/test_had_quadra_batched_directional_scaffold.cpp
  run_had_quadra_batched_directional_scaffold_test.sh

Run:
  ./run_had_quadra_batched_directional_scaffold_test.sh

Expected:
  had_quadra batched directional scaffold tests passed

EOF
