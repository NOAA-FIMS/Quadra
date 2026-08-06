#include "../core/laplace/random_effect_newton.hpp"

#include <cmath>
#include <iostream>
#include <vector>

int main() {
  quadra::ParameterPartition partition;
  partition.fixed_indices_m = {0};
  partition.random_indices_m = {1};
  quadra::RandomEffectNewtonOptions options;
  options.max_iterations_m = 5;
  options.gradient_tolerance_m = 1.0e-10;
  options.step_tolerance_m = 1.0e-12;

  auto evaluate = [](const std::vector<double> &random) {
    quadra::RandomEffectHessianResult out;
    const double u = random[0];
    // Mimic replay/evaluation noise at a large objective: the exact Newton
    // candidate improves the gradient to zero, while its reported objective
    // is one picounit higher and therefore fails a strict Armijo test.
    out.objective_value_m = 470.0 + (std::abs(u) < 5.0e-5 ? 1.0e-12 : 0.0);
    out.gradient_random_m = {u};
    out.gradient_norm_m = std::abs(u);
    out.random_m = random;
    out.full_m = {0.0, u};
    out.gradient_fixed_m = {0.0};
    out.mixed_hessian_m = Eigen::MatrixXd::Zero(1, 1);
    out.hessian_random_m.resize(1, 1);
    out.hessian_random_m.insert(0, 0) = 1.0;
    out.hessian_random_m.makeCompressed();
    return out;
  };

  const auto result = quadra::optimize_random_effects_newton_with_evaluator(
      {0.0}, {1.0e-4}, partition, options, evaluate);
  if (!result.converged_m || result.gradient_norm_m > 1.0e-10 ||
      std::abs(result.u_hat_m[0]) > 1.0e-12) {
    std::cerr << "FAIL: roundoff-scale gradient-improving step rejected\n";
    return 1;
  }
  std::cout << "PASS: Newton accepts gradient improvement within objective "
               "roundoff\n";
  return 0;
}
