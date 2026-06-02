#!/usr/bin/env bash
set -euo pipefail

# install_had_graph_workspace_v1.sh
#
# Adds first-class reusable HAD graph workspace infrastructure:
#
#   core/had_graph_workspace.hpp
#
# This is intentionally small:
#   - owns a had::ADGraph
#   - activates it through had::g_ADGraph
#   - supports Build(builder)
#   - supports ResetAdjoints()
#   - supports PropagateAdjoint(output)
#   - supports ResizeDirectionalBatch(K)
#
# It does not yet change Laplace/exact-gradient production code.

mkdir -p core tests .quadra_patch_backups

backup_file() {
  local f="$1"
  if [[ -f "$f" ]]; then
    local b=".quadra_patch_backups/$(echo "$f" | tr '/ ' '__').$(date +%Y%m%d_%H%M%S).bak"
    cp "$f" "$b"
    echo "Backed up $f -> $b"
  fi
}

backup_file core/had_graph_workspace.hpp
backup_file tests/test_had_graph_workspace.cpp
backup_file run_had_graph_workspace_test.sh

if [[ ! -f core/had_quadra.hpp ]]; then
  echo "ERROR: missing core/had_quadra.hpp"
  exit 1
fi

cat > core/had_graph_workspace.hpp <<'EOF'
#pragma once

#include <cstddef>
#include <stdexcept>
#include <utility>

#include "had_quadra.hpp"

namespace quadra {

// Owns and reuses a HAD graph.
//
// This is the production-facing seam for graph reuse. It intentionally keeps
// the first API small and explicit:
//
//   HadGraphWorkspace ws;
//   auto y = ws.Build([&]() { return model(theta, u); });
//   ws.PropagateAdjoint(y);
//
// Later layers can build exact-gradient/Laplace workspaces on top of this
// without directly managing had::g_ADGraph.
class HadGraphWorkspace {
public:
    HadGraphWorkspace() = default;

    HadGraphWorkspace(const HadGraphWorkspace&) = delete;
    HadGraphWorkspace& operator=(const HadGraphWorkspace&) = delete;

    HadGraphWorkspace(HadGraphWorkspace&&) = default;
    HadGraphWorkspace& operator=(HadGraphWorkspace&&) = default;

    had::ADGraph& Graph() { return graph_; }
    const had::ADGraph& Graph() const { return graph_; }

    void Activate() { had::g_ADGraph = &graph_; }

    bool IsActive() const { return had::g_ADGraph == &graph_; }

    std::size_t VertexCount() const { return graph_.vertices.size(); }

    void Clear() {
        graph_ = had::ADGraph();
        built_ = false;
        output_var_id_ = 0;
    }

    template <class Builder>
    had::AReal Build(Builder&& builder) {
        Clear();
        Activate();

        had::AReal output = std::forward<Builder>(builder)();

        built_ = true;
        output_var_id_ = output.varId;

        return output;
    }

    void ResetAdjoints() {
        Activate();

        for (auto& vertex : graph_.vertices) {
            vertex.w = had::Real(0.0);
            vertex.wDot = had::Real(0.0);
            vertex.soWDot = had::Real(0.0);
            vertex.dot = had::Real(0.0);

            if (!vertex.wDotBatch.empty()) {
                std::fill(vertex.wDotBatch.begin(),
                          vertex.wDotBatch.end(),
                          had::Real(0.0));
            }
            if (!vertex.soWDotBatch.empty()) {
                std::fill(vertex.soWDotBatch.begin(),
                          vertex.soWDotBatch.end(),
                          had::Real(0.0));
            }
            if (!vertex.dotBatch.empty()) {
                std::fill(vertex.dotBatch.begin(),
                          vertex.dotBatch.end(),
                          had::Real(0.0));
            }
            if (!vertex.e1.dwBatch.empty()) {
                std::fill(vertex.e1.dwBatch.begin(),
                          vertex.e1.dwBatch.end(),
                          had::Real(0.0));
            }
            if (!vertex.e2.dwBatch.empty()) {
                std::fill(vertex.e2.dwBatch.begin(),
                          vertex.e2.dwBatch.end(),
                          had::Real(0.0));
            }
        }

        if (graph_.soEdges.size() < graph_.vertices.size()) {
            graph_.soEdges.resize(graph_.vertices.size());
        } else {
            for (auto& tree : graph_.soEdges) {
                tree.Clear();
            }
        }

        if (graph_.soEdgesDot.size() < graph_.vertices.size()) {
            graph_.soEdgesDot.resize(graph_.vertices.size());
        } else {
            for (auto& tree : graph_.soEdgesDot) {
                tree.Clear();
            }
        }

        graph_.selfSoEdges.assign(graph_.vertices.size(), had::Real(0.0));
        graph_.selfSoEdgesDot.assign(graph_.vertices.size(), had::Real(0.0));
    }

    void PropagateAdjoint() {
        if (!built_) {
            throw std::logic_error(
                "HadGraphWorkspace::PropagateAdjoint called before Build.");
        }

        PropagateAdjoint(output_var_id_);
    }

    void PropagateAdjoint(had::VertexId output_var_id) {
        Activate();
        ResetAdjoints();

        if (output_var_id >= graph_.vertices.size()) {
            throw std::out_of_range(
                "HadGraphWorkspace::PropagateAdjoint output_var_id out of range.");
        }

        graph_.vertices[output_var_id].w = had::Real(1.0);
        had::PropagateAdjoint();
    }

    void ResizeDirectionalBatch(std::size_t n_directions) {
        Activate();
        had::ResizeDirectionalBatch(static_cast<int>(n_directions));
    }

    void PropagateAdjointDirectionalBatch() {
        Activate();
        had::PropagateAdjointDirectionalBatch();
    }

    had::VertexId OutputVarId() const {
        if (!built_) {
            throw std::logic_error(
                "HadGraphWorkspace::OutputVarId called before Build.");
        }
        return output_var_id_;
    }

private:
    had::ADGraph graph_;
    bool built_ = false;
    had::VertexId output_var_id_ = 0;
};

}  // namespace quadra
EOF

cat > tests/test_had_graph_workspace.cpp <<'EOF'
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "../core/had_graph_workspace.hpp"

DECLARE_ADGRAPH()

namespace {

had::AReal objective(const had::AReal& x, const had::AReal& y) {
    return exp(x * y) + x * x + sin(y);
}

}  // namespace

int main() {
    quadra::HadGraphWorkspace workspace;

    had::AReal x;
    had::AReal y;
    had::AReal f;

    f = workspace.Build([&]() {
        x = had::AReal(0.4);
        y = had::AReal(-0.7);
        return objective(x, y);
    });

    if (!workspace.IsActive()) {
        throw std::runtime_error("workspace did not activate HAD graph");
    }

    if (workspace.VertexCount() == 0) {
        throw std::runtime_error("workspace graph has no vertices");
    }

    workspace.PropagateAdjoint(f.varId);

    const double hxy_scalar = had::GetAdjoint(x, y);

    workspace.ResizeDirectionalBatch(2);
    had::SetARealDotBatch(x, 0, 1.0);
    had::SetARealDotBatch(y, 0, 0.0);
    had::SetARealDotBatch(x, 1, 0.0);
    had::SetARealDotBatch(y, 1, 1.0);

    workspace.PropagateAdjointDirectionalBatch();

    const double hxy_batch_0 = had::GetAdjointDotBatch(x, y, 0);
    const double hxy_batch_1 = had::GetAdjointDotBatch(x, y, 1);

    // The exact values are not important for this workspace test. We only need
    // to prove that the graph-owned batch path is live and queryable.
    if (!std::isfinite(hxy_scalar) ||
        !std::isfinite(hxy_batch_0) ||
        !std::isfinite(hxy_batch_1)) {
        throw std::runtime_error("workspace produced non-finite derivative values");
    }

    workspace.ResetAdjoints();

    if (workspace.Graph().vertices[f.varId].w != 0.0) {
        throw std::runtime_error("ResetAdjoints did not clear output adjoint");
    }

    std::cout << "had graph workspace tests passed\n";
    std::cout << "vertices = " << workspace.VertexCount() << "\n";
    std::cout << "hxy scalar = " << hxy_scalar << "\n";
    std::cout << "hxy batch 0 = " << hxy_batch_0 << "\n";
    std::cout << "hxy batch 1 = " << hxy_batch_1 << "\n";

    return 0;
}
EOF

cat > run_had_graph_workspace_test.sh <<'EOF'
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
  tests/test_had_graph_workspace.cpp \
  -o build/tests/test_had_graph_workspace

./build/tests/test_had_graph_workspace
EOF

chmod +x run_had_graph_workspace_test.sh

cat <<'EOF'

Installed HadGraphWorkspace v1.

Files added:
  core/had_graph_workspace.hpp
  tests/test_had_graph_workspace.cpp
  run_had_graph_workspace_test.sh

Run:
  ./run_had_graph_workspace_test.sh
  ./run_had_quadra_directional_batch_propagation_test.sh

EOF
