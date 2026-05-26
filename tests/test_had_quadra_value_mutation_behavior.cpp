#include <cmath>
#include <iostream>

#include "../core/autodiff.hpp"

DECLARE_ADGRAPH();

// This is the critical behavioral test for reusable tapes.
//
// Goal:
//
// 1. Build ONE graph/tape.
// 2. Construct f(x) = x^2.
// 3. Run reverse pass.
// 4. Mutate x.val in-place.
// 5. Re-evaluate/reverse WITHOUT rebuilding the graph.
// 6. Verify objective and gradient update correctly.
//
// If this works, Quadra can likely support reusable tapes without rebuilding
// graphs every optimizer iteration.

int main() {

  std::cout << "had_quadra value mutation behavior test\n";

  quadra::TapeContext tape;
  quadra::ADScope scope(tape.graph);

  quadra::AD x = 2.0;

  quadra::AD y = x * x;

  // First reverse pass
  scope.backward(y);

  double value1 = quadra::value_of(y);

  Eigen::VectorXd g1 = quadra::extract_gradient(std::vector<quadra::AD>{x});

  std::cout << "initial value = " << value1 << "\n";
  std::cout << "initial grad = " << g1[0] << "\n";

  if (std::abs(value1 - 4.0) > 1e-12) {
    std::cout << "FAIL: initial value mismatch\n";
    return 1;
  }

  if (std::abs(g1[0] - 4.0) > 1e-12) {
    std::cout << "FAIL: initial gradient mismatch\n";
    return 1;
  }

  // ---------------------------------------------------------------------
  // MUTATE VALUE IN PLACE
  // ---------------------------------------------------------------------

  x.val = 3.0;

  // IMPORTANT:
  //
  // We intentionally DO NOT rebuild:
  //
  //   quadra::AD y = x * x;
  //
  // The test is checking whether the existing graph updates.

  scope.backward(y);

  double value2 = quadra::value_of(y);

  Eigen::VectorXd g2 = quadra::extract_gradient(std::vector<quadra::AD>{x});

  std::cout << "mutated value = " << value2 << "\n";
  std::cout << "mutated grad = " << g2[0] << "\n";

  const bool value_updated = std::abs(value2 - 9.0) <= 1e-12;

  const bool gradient_updated = std::abs(g2[0] - 6.0) <= 1e-12;

  std::cout << "value updated correctly = " << value_updated << "\n";

  std::cout << "gradient updated correctly = " << gradient_updated << "\n";

  if (!value_updated || !gradient_updated) {

    std::cout << "\n";
    std::cout << "RESULT:\n";
    std::cout << "In-place mutation does NOT fully update the graph.\n";
    std::cout << "Reusable tapes will require explicit forward propagation\n";
    std::cout << "or graph-node update support in had_quadra.\n";

    return 1;
  }

  std::cout << "\n";
  std::cout << "SUCCESS:\n";
  std::cout << "Existing graph responds to in-place AD value mutation.\n";
  std::cout << "Reusable tape evaluators are feasible.\n";

  return 0;
}
