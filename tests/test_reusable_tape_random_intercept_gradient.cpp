#include <cmath>
#include <iostream>
#include <vector>

#include "../core/autodiff.hpp"

DECLARE_ADGRAPH();

template <typename Type>
Type random_intercept_joint(const std::vector<double> &y, const Type &mu,
                            const Type &u) {
  Type nll = Type(0.0);

  for (double yi : y) {
    Type r = Type(yi) - (mu + u);
    nll += Type(0.5) * r * r;
  }

  nll += Type(0.5) * u * u;

  return nll;
}

int main() {
  std::cout << "Reusable tape random-intercept gradient test\n";

  std::vector<double> y = {4.8, 5.1, 5.0, 4.9, 5.2};

  quadra::TapeContext tape;
  quadra::ADScope scope(tape.graph);

  quadra::AD mu = 4.7;
  quadra::AD u = 0.25;

  quadra::AD nll = random_intercept_joint(y, mu, u);

  scope.backward(nll);

  Eigen::VectorXd g1 = quadra::extract_gradient(std::vector<quadra::AD>{mu, u});

  double f1 = quadra::value_of(nll);

  std::cout << "initial f = " << f1 << "\n";
  std::cout << "initial grad_mu = " << g1[0] << "\n";
  std::cout << "initial grad_u = " << g1[1] << "\n";

  // Move to the optimum for this toy model.
  quadra::set_value(mu, 5.0);
  quadra::set_value(u, 0.0);

  scope.forward();
  scope.zero_adjoints();
  scope.backward(nll);

  Eigen::VectorXd g2 = quadra::extract_gradient(std::vector<quadra::AD>{mu, u});

  double f2 = quadra::value_of(nll);

  std::cout << "replayed f = " << f2 << "\n";
  std::cout << "replayed grad_mu = " << g2[0] << "\n";
  std::cout << "replayed grad_u = " << g2[1] << "\n";

  // At mu=5, u=0:
  // sum(y)=25, n=5
  // grad_mu = n*(mu+u)-sum(y)=0
  // grad_u = n*(mu+u)-sum(y)+u=0
  if (std::abs(g2[0]) > 1e-10) {
    std::cerr << "FAIL: replayed grad_mu mismatch\n";
    return 1;
  }

  if (std::abs(g2[1]) > 1e-10) {
    std::cerr << "FAIL: replayed grad_u mismatch\n";
    return 1;
  }

  std::cout << "PASS\n";
  return 0;
}
