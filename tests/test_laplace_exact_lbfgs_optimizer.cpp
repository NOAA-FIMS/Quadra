#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "../core/laplace/laplace_exact_lbfgs_optimizer.hpp"
#include "../core/model/quadra_model.hpp"

DECLARE_ADGRAPH();

class RandomInterceptModel : public quadra::QuadraModel<RandomInterceptModel> {
public:
  RandomInterceptModel(std::vector<double> y) : y_m(std::move(y)) {
    parameters_m.add("mu", 0.0, quadra::ParameterTransform::Identity, false);
    parameters_m.add("u", 0.0, quadra::ParameterTransform::Identity, true);
  }

  std::vector<std::string> parameter_names_impl() const {
    return parameters_m.names();
  }

  const quadra::ParameterSet &parameters() const { return parameters_m; }

  template <typename Type>
  Type evaluate_impl(const std::vector<Type> &p,
                     quadra::ModelReportContext &) const {
    Type mu = p[0];
    Type u = p[1];
    Type nll = Type(0.0);

    for (double yi : y_m) {
      Type r = Type(yi) - (mu + u);
      nll += Type(0.5) * r * r;
    }

    nll += Type(0.5) * u * u;
    return nll;
  }

private:
  std::vector<double> y_m;
  quadra::ParameterSet parameters_m;
};

int main() {
  std::cout << "Testing exact-gradient Laplace LBFGS optimizer\n";

  RandomInterceptModel model({4.8, 5.1, 5.0, 4.9, 5.2});

  quadra::LaplaceExactLBFGSOptions options;
  options.max_iterations_m = 100;
  options.memory_m = 5;
  options.gradient_tolerance_m = 1e-8;
  options.step_tolerance_m = 1e-12;
  options.initial_step_scale_m = 1.0;
  options.gradient_m.objective_m.include_constant_m = true;
  options.gradient_m.objective_m.newton_m.gradient_tolerance_m = 1e-12;
  options.gradient_m.objective_m.newton_m.step_tolerance_m = 1e-14;

  auto result = quadra::optimize_laplace_fixed_effects_exact_lbfgs(
      model, {4.0}, {0.0}, model.parameters(), options);

  std::cout << "converged = " << result.converged_m << "\n";
  std::cout << "message = " << result.message_m << "\n";
  std::cout << "iterations = " << result.iterations_m << "\n";
  std::cout << "theta_hat = " << result.theta_hat_m[0] << "\n";
  std::cout << "u_hat = " << result.u_hat_m[0] << "\n";
  std::cout << "laplace objective = " << result.laplace_objective_m << "\n";
  std::cout << "gradient norm = " << result.gradient_norm_m << "\n";

  if (!result.converged_m)
    return 1;
  if (std::abs(result.theta_hat_m[0] - 5.0) > 1e-6)
    return 1;
  if (std::abs(result.u_hat_m[0]) > 1e-6)
    return 1;
  if (result.gradient_norm_m > 1e-7)
    return 1;

  std::cout << "PASS\n";
  return 0;
}
