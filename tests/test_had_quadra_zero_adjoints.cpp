#include <cmath>
#include <iostream>
#include <vector>

#include "../core/autodiff.hpp"

DECLARE_ADGRAPH();

int main() {
  std::cout << "had_quadra zero_adjoints behavior test\n";

  quadra::TapeContext tape;
  quadra::ADScope scope(tape.graph);

  quadra::AD x = 2.0;
  quadra::AD y = x * x;

  scope.backward(y);
  Eigen::VectorXd g1 = quadra::extract_gradient(std::vector<quadra::AD>{x});

  scope.zero_adjoints();

  scope.backward(y);
  Eigen::VectorXd g2 = quadra::extract_gradient(std::vector<quadra::AD>{x});

  std::cout << "grad after first backward = " << g1[0] << "\n";
  std::cout << "grad after zero + second backward = " << g2[0] << "\n";

  if (std::abs(g1[0] - 4.0) > 1e-12) {
    std::cerr << "FAIL: first gradient mismatch\n";
    return 1;
  }

  if (std::abs(g2[0] - 4.0) > 1e-12) {
    std::cerr << "FAIL: zero_adjoints did not prevent accumulation\n";
    return 1;
  }

  quadra::zero_adjoints(tape);

  scope.backward(y);
  Eigen::VectorXd g3 = quadra::extract_gradient(std::vector<quadra::AD>{x});

  std::cout << "grad after free-wrapper zero + third backward = " << g3[0]
            << "\n";

  if (std::abs(g3[0] - 4.0) > 1e-12) {
    std::cerr
        << "FAIL: free-wrapper zero_adjoints did not prevent accumulation\n";
    return 1;
  }

  std::cout << "PASS\n";
  return 0;
}
