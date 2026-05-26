#include <cmath>
#include <iostream>
#include <vector>

#include "../core/autodiff.hpp"

DECLARE_ADGRAPH();

double f(double x) {
  return std::exp(x) + std::log(x) + std::sqrt(x) + 10.0 / x + x / 2.0;
}

double df(double x) {
  return std::exp(x) + 1.0 / x + 0.5 / std::sqrt(x) - 10.0 / (x * x) + 0.5;
}

int main() {
  std::cout << "had_quadra forward replay math-ops test\n";

  quadra::TapeContext tape;
  quadra::ADScope scope(tape.graph);

  quadra::AD x = 2.0;

  quadra::AD y = exp(x) + log(x) + sqrt(x) + 10.0 / x + x / 2.0;

  scope.backward(y);

  double y1 = quadra::value_of(y);

  Eigen::VectorXd g1 = quadra::extract_gradient(std::vector<quadra::AD>{x});

  std::cout << "initial y = " << y1 << "\n";
  std::cout << "initial grad = " << g1[0] << "\n";

  if (std::abs(y1 - f(2.0)) > 1e-10) {
    std::cerr << "FAIL: initial value mismatch\n";
    return 1;
  }

  if (std::abs(g1[0] - df(2.0)) > 1e-10) {
    std::cerr << "FAIL: initial gradient mismatch\n";
    return 1;
  }

  quadra::set_value(x, 3.0);

  scope.forward();
  scope.zero_adjoints();
  scope.backward(y);

  double y2 = quadra::value_of(y);

  Eigen::VectorXd g2 = quadra::extract_gradient(std::vector<quadra::AD>{x});

  std::cout << "replayed y = " << y2 << "\n";
  std::cout << "replayed grad = " << g2[0] << "\n";

  if (std::abs(y2 - f(3.0)) > 1e-10) {
    std::cerr << "FAIL: replayed value mismatch\n";
    std::cerr << "expected " << f(3.0) << "\n";
    return 1;
  }

  if (std::abs(g2[0] - df(3.0)) > 1e-10) {
    std::cerr << "FAIL: replayed gradient mismatch\n";
    std::cerr << "expected " << df(3.0) << "\n";
    return 1;
  }

  std::cout << "PASS\n";
  return 0;
}
