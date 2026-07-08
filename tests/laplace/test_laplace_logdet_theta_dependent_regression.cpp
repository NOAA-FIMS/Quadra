#include "../../core/optimizer.hpp"
#include "../../model/parameter.hpp"

#include <cmath>
#include <iostream>

threadDefine had::ADGraph *had::g_ADGraph = nullptr;

using namespace quadra;

struct ThetaDependentGaussianREModel {
  template <class T>
  T operator()(const std::vector<T> &p) const {
    const T theta = p[0];
    const T u = p[1];

    const T target = T(1.5);
    const T precision = exp(theta);

    return T(0.5) * (theta - target) * (theta - target) +
           T(0.5) * precision * u * u;
  }
};

int main() {
  ParameterVector params;
  params.add(Parameter("theta", 1.25, ParameterTransform::Identity, false));
  params.add(Parameter("u", 0.0, ParameterTransform::Identity, true));

  ThetaDependentGaussianREModel model;

  LaplaceOptions options = default_laplace_options();  options.hessian_drop_tol = 0.0;

  auto result = optimize_lbfgs(model, params, options);

  const double theta = result.par[0];
  const double got = result.fixed_gradient[0];

  // Profiled objective:
  //   0.5(theta - 1.5)^2 + 0.5 log(exp(theta)) - const
  // = 0.5(theta - 1.5)^2 + 0.5 theta - const
  //
  // Gradient:
  //   theta - 1.5 + 0.5
  //
  // Optimum:
  //   theta = 1.0
  const double expected_theta = 1.0;
  const double expected_grad = theta - 1.5 + 0.5;

  const double theta_diff = std::abs(theta - expected_theta);
  const double grad_diff = std::abs(got - expected_grad);

  std::cout << "theta," << theta << "\n";
  std::cout << "theta_expected," << expected_theta << "\n";
  std::cout << "theta_abs_diff," << theta_diff << "\n";
  std::cout << "gradient_got," << got << "\n";
  std::cout << "gradient_expected," << expected_grad << "\n";
  std::cout << "gradient_abs_diff," << grad_diff << "\n";

  if (theta_diff > 1.0e-5 || grad_diff > 1.0e-5) {
    std::cerr << "FAILED: theta-dependent Laplace regression\n";
    return 1;
  }

  std::cout << "PASSED: theta-dependent Laplace regression\n";
  return 0;
}
