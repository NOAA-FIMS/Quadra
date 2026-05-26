#include "../examples/big/catch_at_age_shared.hpp"

#include <cmath>
#include <iostream>

int main() {
  example::CatchAtAgeLaplaceModel model;
  auto params = example::make_big_laplace_parameter_vector();

  quadra::LaplaceOptions opts = quadra::default_laplace_options();

  auto fit = quadra::optimize_lbfgs(model, params, opts);

  const bool ok = std::isfinite(fit.value) && fit.value < 0.0 &&
                  fit.value > -500.0 && fit.grad_norm < 1.0e-3;

  if (!ok) {
    std::cerr << "FAIL: big Laplace black-box convergence contract failed\n";
    std::cerr << "  fit value: " << fit.value << "\n";
    std::cerr << "  fixed gradient norm: " << fit.grad_norm << "\n";
    return 1;
  }

  std::cout << "PASS: big Laplace black-box convergence contract satisfied\n";
  std::cout << "  fit value: " << fit.value << "\n";
  std::cout << "  fixed gradient norm: " << fit.grad_norm << "\n";
  return 0;
}
