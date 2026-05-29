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
