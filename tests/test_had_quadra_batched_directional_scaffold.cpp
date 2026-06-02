#include <iostream>
#include <stdexcept>

#include "../core/had_quadra.hpp"

DECLARE_ADGRAPH()

int main() {
  had::ADGraph graph;
  had::g_ADGraph = &graph;

  had::AReal x(2.0);
  had::AReal y(3.0);

  had::ResizeDirectionalBatch(3);

  had::SetARealDotBatch(x, 0, 1.0);
  had::SetARealDotBatch(x, 1, 2.0);
  had::SetARealDotBatch(y, 2, 4.0);

  if (had::GetVertexDotBatch(x.varId, 0) != 1.0) {
    throw std::runtime_error("x direction 0 was not stored");
  }

  if (had::GetVertexDotBatch(x.varId, 1) != 2.0) {
    throw std::runtime_error("x direction 1 was not stored");
  }

  if (had::GetVertexDotBatch(y.varId, 2) != 4.0) {
    throw std::runtime_error("y direction 2 was not stored");
  }

  had::ClearDirectionalBatch();

  if (had::GetVertexDotBatch(x.varId, 0) != 0.0) {
    throw std::runtime_error("clear did not reset x direction 0");
  }

  if (had::GetAdjointDotBatch(x, y, 0) != 0.0) {
    throw std::runtime_error("empty batched adjoint dot should be zero");
  }

  bool threw = false;
  try {
    had::SetARealDotBatch(x, 4, 1.0);
  } catch (const std::out_of_range &) {
    threw = true;
  }

  if (!threw) {
    throw std::runtime_error("out-of-range batch index did not throw");
  }

  std::cout << "had_quadra batched directional scaffold tests passed\n";
  return 0;
}
