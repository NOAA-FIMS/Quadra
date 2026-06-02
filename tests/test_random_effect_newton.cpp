#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "../core/laplace/random_effect_newton.hpp"
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

    // Observation model: y_i ~ N(mu + u, 1), constants omitted.
    for (double yi : y_m) {
      Type r = Type(yi) - (mu + u);
      nll += Type(0.5) * r * r;
    }

    // Random-effect prior: u ~ N(0, 1), constants omitted.
    nll += Type(0.5) * u * u;

    ctx.report("mu", mu);
    ctx.report("u", u);
    ctx.adreport("mean_prediction", mu + u);

    return nll;
  }

private:
  std::vector<double> y_m;
  quadra::ParameterSet parameters_m;
};

int main() {
  std::cout << "Testing random-effect Newton optimizer\n";

  RandomInterceptModel model({4.8, 5.1, 5.0, 4.9, 5.2});

  const auto partition = quadra::partition_parameters(model.parameters());

  std::vector<double> fixed = {4.7};
  std::vector<double> random_initial = {0.0};

  quadra::RandomEffectNewtonOptions options;
  options.max_iterations_m = 10;
  options.gradient_tolerance_m = 1e-10;
  options.step_tolerance_m = 1e-12;

  auto result = quadra::optimize_random_effects_newton(
      model, fixed, random_initial, partition, options);

  // Analytic solution:
  //
  // gradient = sum(mu + u - y_i) + u
  //          = n * (mu + u) - sum(y) + u
  //
  // Set to zero:
  //
  //   (n + 1)u = sum(y) - n*mu
  //
  // Here sum(y) = 25.0, n = 5, mu = 4.7:
  //
  //   u_hat = (25 - 23.5) / 6 = 0.25
  const double expected_u_hat = 0.25;
  const double expected_hessian = 6.0;

  Eigen::MatrixXd H(result.hessian_random_m);

  std::cout << "converged = " << result.converged_m << "\n";
  std::cout << "message = " << result.message_m << "\n";
  std::cout << "iterations = " << result.iterations_m << "\n";
  std::cout << "u_hat = " << result.u_hat_m[0] << "\n";
  std::cout << "gradient norm = " << result.gradient_norm_m << "\n";
  std::cout << "H_uu = " << H(0, 0) << "\n";

  if (!result.converged_m) {
    std::cerr << "FAIL: Newton optimizer did not converge\n";
    return 1;
  }

  if (std::abs(result.u_hat_m[0] - expected_u_hat) > 1e-10) {
    std::cerr << "FAIL: u_hat mismatch\n";
    return 1;
  }

  if (result.gradient_norm_m > 1e-10) {
    std::cerr << "FAIL: gradient norm too large\n";
    return 1;
  }

  if (std::abs(H(0, 0) - expected_hessian) > 1e-12) {
    std::cerr << "FAIL: Hessian mismatch\n";
    return 1;
  }

  if (result.full_m.size() != 2 || std::abs(result.full_m[0] - 4.7) > 1e-12 ||
      std::abs(result.full_m[1] - expected_u_hat) > 1e-10) {
    std::cerr << "FAIL: full merged vector mismatch\n";
    std::cerr << "full[0] = " << result.full_m[0]
              << ", full[1] = " << result.full_m[1] << "\n";
    return 1;
  }

  std::cout << "PASS\n";
  return 0;
}
