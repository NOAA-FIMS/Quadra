#include "../core/laplace.hpp"

#include <algorithm>
#include <iostream>
#include <utility>
#include <vector>

DECLARE_ADGRAPH();

namespace {

quadra::SparseHessianPattern discover(bool coupled, bool zero_at_point) {
  quadra::TapeContext tape;
  quadra::ADScope scope(tape.graph);
  std::vector<quadra::AD> x = quadra::to_ad(std::vector<double>{0.0, 0.0});

  quadra::AD objective = x[0] * x[0] + x[1] * x[1];
  if (coupled) {
    objective = objective + (zero_at_point ? x[0] * x[0] * x[1] : x[0] * x[1]);
  }
  scope.backward(objective);
  return quadra::get_pattern(scope, x, {0, 1});
}

bool contains(const quadra::SparseHessianPattern &pattern, int row, int col) {
  return std::find(pattern.begin(), pattern.end(), std::make_pair(row, col)) !=
         pattern.end();
}

} // namespace

int main() {
  const auto diagonal = discover(false, false);
  const auto coupled = discover(true, false);
  const auto zero_at_point = discover(true, true);

  if (contains(diagonal, 0, 1) || !contains(coupled, 0, 1)) {
    std::cerr << "same-sized models incorrectly shared a Hessian pattern\n";
    return 1;
  }
  if (!contains(zero_at_point, 0, 1) || !contains(zero_at_point, 1, 0)) {
    std::cerr << "zero-valued structural Hessian edge was dropped\n";
    return 1;
  }

  std::cout << "PASS: Hessian structure discovery is model-safe and "
               "value-independent\n";
  return 0;
}
