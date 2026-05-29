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

    had::g_ADGraph->vertices[f.varId].w = 1.0;
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

    had::g_ADGraph->vertices[f.varId].w = 1.0;
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
