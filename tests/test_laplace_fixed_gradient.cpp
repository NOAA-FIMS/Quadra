#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "../core/laplace/laplace_fixed_gradient.hpp"
#include "../core/model/parameter_partition.hpp"
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
                     quadra::ModelReportContext &ctx) const {
    Type mu = p[0];
    Type u = p[1];

    Type nll = Type(0.0);

    for (double yi : y_m) {
      Type r = Type(yi) - (mu + u);
      nll += Type(0.5) * r * r;
    }

    nll += Type(0.5) * u * u;

    ctx.report("mu", mu);
    ctx.report("u_hat", u);
    ctx.adreport("mean_prediction", mu + u);

    return nll;
  }

private:
  std::vector<double> y_m;
  quadra::ParameterSet parameters_m;
};

int main() {
  std::cout << "Testing Laplace fixed-effect gradient\n";

  RandomInterceptModel model({4.8, 5.1, 5.0, 4.9, 5.2});

  std::vector<double> fixed = {4.7};
  std::vector<double> u0 = {0.0};

  quadra::LaplaceFixedGradientOptions options;
  options.relative_step_m = 1e-6;
  options.absolute_step_m = 1e-7;
  options.use_central_difference_m = true;
  options.objective_m.include_constant_m = true;
  options.objective_m.newton_m.gradient_tolerance_m = 1e-12;
  options.objective_m.newton_m.step_tolerance_m = 1e-14;

  auto result = quadra::evaluate_laplace_fixed_gradient(
      model, fixed, u0, model.parameters(), options);

  // Analytic marginal derivative for this Gaussian random-intercept model:
  //
  // f(mu, u) = 0.5 * sum_i (y_i - (mu + u))^2 + 0.5 * u^2
  //
  // u_hat(mu) = (sum(y) - n*mu) / (n + 1)
  //
  // By the envelope theorem, d/dmu f(mu, u_hat(mu)) =
  // partial f / partial mu at u_hat =
  // sum_i (mu + u_hat - y_i).
  //
  // At mu = 4.7, u_hat = 0.25, mu + u_hat = 4.95,
  // sum_i (4.95 - y_i) = 5 * 4.95 - 25 = -0.25.
  //
  // logdet(H_uu) is constant because H_uu = n + 1, so its derivative is 0.
  const double expected_gradient = -0.25;
  const double expected_u_hat = 0.25;

  std::cout << "converged = " << result.converged_m << "\n";
  std::cout << "logdet ok = " << result.logdet_ok_m << "\n";
  std::cout << "laplace objective = " << result.laplace_objective_m << "\n";
  std::cout << "u_hat = " << result.u_hat_m[0] << "\n";
  std::cout << "gradient_mu = " << result.gradient_fixed_m[0] << "\n";
  std::cout << "gradient norm = " << result.gradient_norm_m << "\n";

  if (!result.converged_m || !result.logdet_ok_m) {
    std::cerr << "FAIL: base Laplace objective failed\n";
    return 1;
  }

  if (std::abs(result.u_hat_m[0] - expected_u_hat) > 1e-10) {
    std::cerr << "FAIL: u_hat mismatch\n";
    return 1;
  }

  if (std::abs(result.gradient_fixed_m[0] - expected_gradient) > 1e-5) {
    std::cerr << "FAIL: fixed-effect gradient mismatch\n";
    return 1;
  }

  std::cout << "PASS\n";
  return 0;
}
