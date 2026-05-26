#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "../core/autodiff/model_gradient.hpp"
#include "../core/model/parameter.hpp"
#include "../core/model/quadra_model.hpp"

DECLARE_ADGRAPH();

class TransformedNormalModel
    : public quadra::QuadraModel<TransformedNormalModel> {
public:
  TransformedNormalModel() {
    parameters_m.add("log_sigma", 0.0, quadra::ParameterTransform::Log);
    parameters_m.add("logit_p", 0.0, quadra::ParameterTransform::Logit);
    parameters_m.add("squared_q", 2.0, quadra::ParameterTransform::Square);
  }

  std::vector<std::string> parameter_names_impl() const {
    return parameters_m.names();
  }

  std::vector<quadra::ParameterTransform> parameter_transforms() const {
    return parameters_m.transforms();
  }

  template <typename Type>
  Type evaluate_impl(const std::vector<Type> &unconstrained,
                     quadra::ModelReportContext &ctx) const {
    std::vector<Type> p =
        quadra::apply_transforms(unconstrained, parameters_m.transforms());

    Type sigma = p[0];
    Type prob = p[1];
    Type q = p[2];

    ctx.report("sigma", sigma);
    ctx.report("probability", prob);
    ctx.report("q", q);

    // Simple objective with minimum at sigma = 1, prob = 0.5, q = 4.
    Type nll = Type(0.0);
    nll += Type(0.5) * (sigma - Type(1.0)) * (sigma - Type(1.0));
    nll += Type(0.5) * (prob - Type(0.5)) * (prob - Type(0.5));
    nll += Type(0.5) * (q - Type(4.0)) * (q - Type(4.0));

    return nll;
  }

private:
  quadra::ParameterSet parameters_m;
};

int main() {
  std::cout << "Testing parameter transforms\n";

  TransformedNormalModel model;
  std::vector<double> x = {0.0, 0.0, 2.0};

  auto result = quadra::evaluate_gradient(model, x);

  std::cout << "objective = " << result.objective_value_m << "\n";
  std::cout << "gradient norm = " << result.gradient_norm_m << "\n";

  for (const auto &r : result.reports_m) {
    std::cout << r.name_m << " = " << r.value_m << "\n";
  }

  if (std::abs(result.objective_value_m) > 1e-12) {
    std::cerr << "FAIL: expected transformed optimum to have objective zero\n";
    return 1;
  }

  if (std::abs(result.gradient_norm_m) > 1e-12) {
    std::cerr
        << "FAIL: expected transformed optimum to have gradient norm zero\n";
    return 1;
  }

  if (result.reports_m.size() != 3) {
    std::cerr << "FAIL: expected three transformed reports\n";
    return 1;
  }

  const double sigma = result.reports_m[0].value_m;
  const double probability = result.reports_m[1].value_m;
  const double q = result.reports_m[2].value_m;

  if (std::abs(sigma - 1.0) > 1e-12 || std::abs(probability - 0.5) > 1e-12 ||
      std::abs(q - 4.0) > 1e-12) {
    std::cerr << "FAIL: transformed values were incorrect\n";
    return 1;
  }

  std::cout << "PASS\n";
  return 0;
}
