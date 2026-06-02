#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "../core/had_quadra.hpp"

DECLARE_ADGRAPH()

namespace {

void require_finite(const char *name, double x) {
  if (!std::isfinite(x)) {
    std::cerr << name << " is non-finite: " << x << "\n";
    std::exit(1);
  }
}

void require_nonzero(const char *name, double x, double tol = 1.0e-14) {
  require_finite(name, x);
  if (std::abs(x) <= tol) {
    std::cerr << name << " unexpectedly zero: " << x << "\n";
    std::exit(1);
  }
}

void print_vertex(const char *label, const had::AReal &x) {
  const auto &v = had::g_ADGraph->vertices[x.varId];

  std::cout << label << " varId=" << x.varId << " val=" << x.val
            << " dot=" << x.dot << " vertex.dot=" << v.dot;

  if (!x.dotBatch.empty()) {
    std::cout << " areal.dotBatch[0]=" << x.dotBatch[0];
  }
  if (!v.dotBatch.empty()) {
    std::cout << " vertex.dotBatch[0]=" << v.dotBatch[0];
  }
  if (!v.wDotBatch.empty()) {
    std::cout << " wDotBatch[0]=" << v.wDotBatch[0];
  }
  if (!v.e1.dwBatch.empty()) {
    std::cout << " e1.dwBatch[0]=" << v.e1.dwBatch[0];
  }
  if (!v.e2.dwBatch.empty()) {
    std::cout << " e2.dwBatch[0]=" << v.e2.dwBatch[0];
  }

  std::cout << "\n";
}

} // namespace

int main() {
  had::ADGraph graph;
  had::g_ADGraph = &graph;

  had::AReal x(0.7);
  had::AReal y(-0.4);

  // Deliberately nonlinear and mixed:
  //
  //   f = exp(x*y) + x*x*y + sin(y)
  //
  // This should produce non-zero Hessian directional terms for many seeds.
  had::AReal f = exp(x * y) + x * x * y + sin(y);

  had::g_ADGraph->vertices[f.varId].w = 1.0;
  had::PropagateAdjoint();

  had::ResizeDirectionalBatch(2);

  // Direction 0: dx = 1, dy = 0
  had::SetARealDotBatch(x, 0, 1.0);
  had::SetARealDotBatch(y, 0, 0.0);

  // Direction 1: dx = 0, dy = 1
  had::SetARealDotBatch(x, 1, 0.0);
  had::SetARealDotBatch(y, 1, 1.0);

  std::cout << "After seeding batch directions:\n";
  print_vertex("x", x);
  print_vertex("y", y);
  print_vertex("f", f);

  if (x.dotBatch.empty() || y.dotBatch.empty()) {
    throw std::runtime_error("AReal dotBatch was not initialized");
  }

  if (had::g_ADGraph->vertices[x.varId].dotBatch.empty() ||
      had::g_ADGraph->vertices[y.varId].dotBatch.empty()) {
    throw std::runtime_error("vertex dotBatch was not initialized");
  }

  require_nonzero("x.dotBatch[0]", x.dotBatch[0]);
  require_nonzero("vertex x.dotBatch[0]",
                  had::g_ADGraph->vertices[x.varId].dotBatch[0]);
  require_nonzero("y.dotBatch[1]", y.dotBatch[1]);
  require_nonzero("vertex y.dotBatch[1]",
                  had::g_ADGraph->vertices[y.varId].dotBatch[1]);

  // PropagateAdjoint() consumes/clears reverse adjoints while building the
  // base Hessian, so the directional reverse pass needs the output adjoint
  // seed restored.
  had::g_ADGraph->vertices[f.varId].w = 1.0;
  had::PropagateAdjointDirectionalBatch();

  std::cout << "\nAfter PropagateAdjointDirectionalBatch:\n";
  print_vertex("x", x);
  print_vertex("y", y);
  print_vertex("f", f);

  const double hxx_d0 = had::GetAdjointDotBatch(x, x, 0);
  const double hxy_d0 = had::GetAdjointDotBatch(x, y, 0);
  const double hyy_d0 = had::GetAdjointDotBatch(y, y, 0);

  const double hxx_d1 = had::GetAdjointDotBatch(x, x, 1);
  const double hxy_d1 = had::GetAdjointDotBatch(x, y, 1);
  const double hyy_d1 = had::GetAdjointDotBatch(y, y, 1);

  std::cout << "\nBatch Hdot readout:\n";
  std::cout << "direction 0: " << "hxx_dot=" << hxx_d0 << " hxy_dot=" << hxy_d0
            << " hyy_dot=" << hyy_d0 << "\n";
  std::cout << "direction 1: " << "hxx_dot=" << hxx_d1 << " hxy_dot=" << hxy_d1
            << " hyy_dot=" << hyy_d1 << "\n";

  require_finite("hxx_d0", hxx_d0);
  require_finite("hxy_d0", hxy_d0);
  require_finite("hyy_d0", hyy_d0);
  require_finite("hxx_d1", hxx_d1);
  require_finite("hxy_d1", hxy_d1);
  require_finite("hyy_d1", hyy_d1);

  const double signal = std::abs(hxx_d0) + std::abs(hxy_d0) + std::abs(hyy_d0) +
                        std::abs(hxx_d1) + std::abs(hxy_d1) + std::abs(hyy_d1);

  if (signal <= 1.0e-14) {
    std::cerr << "\nNo non-zero batch Hdot signal detected.\n";
    std::cerr << "This means seeding works, but directional batch "
                 "information is not reaching Hdot readout.\n";
    return 2;
  }

  std::cout << "\nhad_quadra nonzero batch directional diagnostic passed\n";
  return 0;
}
