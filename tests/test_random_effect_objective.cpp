#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "../core/laplace/random_effect_objective.hpp"
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
  std::cout << "Testing random-effect objective wrapper\n";

  RandomInterceptModel model({4.8, 5.1, 5.0, 4.9, 5.2});

  const auto partition = quadra::partition_parameters(model.parameters());

  std::vector<double> fixed = {5.0};
  std::vector<double> random = {0.0};

  auto result = quadra::evaluate_random_effect_objective_gradient(
      model, fixed, random, partition);

  // For y centered at 5 and mu = 5, u = 0 is the optimum.
  // Gradient wrt u = sum(mu + u - y_i) + u = 0.
  const double expected_gradient = 0.0;
  const double expected_objective = 0.05;

  std::cout << "objective = " << result.objective_value_m << "\n";
  std::cout << "gradient_u = " << result.gradient_random_m[0] << "\n";
  std::cout << "gradient norm = " << result.gradient_norm_m << "\n";

  if (std::abs(result.objective_value_m - expected_objective) > 1e-12) {
    std::cerr << "FAIL: objective mismatch\n";
    return 1;
  }

  if (std::abs(result.gradient_random_m[0] - expected_gradient) > 1e-12) {
    std::cerr << "FAIL: random-effect gradient mismatch\n";
    return 1;
  }

  if (result.full_m != std::vector<double>{5.0, 0.0}) {
    std::cerr << "FAIL: merged full vector mismatch\n";
    return 1;
  }

  if (result.reports_m.size() != 3) {
    std::cerr << "FAIL: expected three reports\n";
    return 1;
  }

  // Move u away from optimum and check analytic gradient.
  random = {0.2};

  result = quadra::evaluate_random_effect_objective_gradient(model, fixed,
                                                             random, partition);

  // gradient = sum(mu + u - y_i) + u
  //          = 5 * 0.2 + 0.2 = 1.2
  const double expected_gradient_shifted = 1.2;

  std::cout << "shifted gradient_u = " << result.gradient_random_m[0] << "\n";

  if (std::abs(result.gradient_random_m[0] - expected_gradient_shifted) >
      1e-12) {
    std::cerr << "FAIL: shifted random-effect gradient mismatch\n";
    return 1;
  }

  std::cout << "PASS\n";
  return 0;
}
