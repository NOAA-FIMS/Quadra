#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "../core/laplace/random_effect_hessian.hpp"
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
  std::cout << "Testing random-effect sparse Hessian wrapper\n";

  RandomInterceptModel model({4.8, 5.1, 5.0, 4.9, 5.2});

  const auto partition = quadra::partition_parameters(model.parameters());

  std::vector<double> fixed = {5.0};
  std::vector<double> random = {0.0};

  auto result =
      quadra::evaluate_random_effect_hessian(model, fixed, random, partition);

  const double expected_objective = 0.05;
  const double expected_gradient = 0.0;

  // For f(u) = 0.5 * sum_i (y_i - (mu + u))^2 + 0.5 * u^2,
  // d2f / du2 = n + 1. Here n = 5, so H_uu = 6.
  const double expected_hessian = 6.0;

  Eigen::MatrixXd H = quadra::dense_random_hessian(result);

  std::cout << "objective = " << result.objective_value_m << "\n";
  std::cout << "gradient_u = " << result.gradient_random_m[0] << "\n";
  std::cout << "H_uu = " << H(0, 0) << "\n";
  std::cout << "nnz = " << result.hessian_random_m.nonZeros() << "\n";

  if (std::abs(result.objective_value_m - expected_objective) > 1e-12) {
    std::cerr << "FAIL: objective mismatch\n";
    return 1;
  }

  if (std::abs(result.gradient_random_m[0] - expected_gradient) > 1e-12) {
    std::cerr << "FAIL: gradient mismatch\n";
    return 1;
  }

  if (H.rows() != 1 || H.cols() != 1) {
    std::cerr << "FAIL: Hessian dimensions mismatch\n";
    return 1;
  }

  if (std::abs(H(0, 0) - expected_hessian) > 1e-12) {
    std::cerr << "FAIL: Hessian value mismatch\n";
    return 1;
  }

  if (result.full_m != std::vector<double>{5.0, 0.0}) {
    std::cerr << "FAIL: merged full vector mismatch\n";
    return 1;
  }

  // Shift u away from optimum. Hessian should remain 6 for this model.
  random = {0.2};

  result =
      quadra::evaluate_random_effect_hessian(model, fixed, random, partition);

  H = quadra::dense_random_hessian(result);

  const double expected_shifted_gradient = 1.2;

  std::cout << "shifted gradient_u = " << result.gradient_random_m[0] << "\n";
  std::cout << "shifted H_uu = " << H(0, 0) << "\n";

  if (std::abs(result.gradient_random_m[0] - expected_shifted_gradient) >
      1e-12) {
    std::cerr << "FAIL: shifted gradient mismatch\n";
    return 1;
  }

  if (std::abs(H(0, 0) - expected_hessian) > 1e-12) {
    std::cerr << "FAIL: shifted Hessian value mismatch\n";
    return 1;
  }

  std::cout << "PASS\n";
  return 0;
}
