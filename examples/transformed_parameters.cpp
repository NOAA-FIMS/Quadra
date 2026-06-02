#include <iostream>
#include <string>
#include <vector>

#include "../core/autodiff/model_gradient.hpp"
#include "../core/model/parameter.hpp"
#include "../core/model/quadra_model.hpp"

DECLARE_ADGRAPH();

class SimpleTransformedModel
    : public quadra::QuadraModel<SimpleTransformedModel> {
public:
  SimpleTransformedModel() {
    parameters_m.add("log_sigma", -0.25, quadra::ParameterTransform::Log);
    parameters_m.add("logit_p", 0.75, quadra::ParameterTransform::Logit);
    parameters_m.add("sqrt_q", 1.5, quadra::ParameterTransform::Square);
  }

  std::vector<std::string> parameter_names_impl() const {
    return parameters_m.names();
  }

  std::vector<double> initials() const { return parameters_m.initials(); }

  template <typename Type>
  Type evaluate_impl(const std::vector<Type> &unconstrained,
                     quadra::ModelReportContext &ctx) const {
    std::vector<Type> p =
        quadra::apply_transforms(unconstrained, parameters_m.transforms());

    Type sigma = p[0];
    Type p_survive = p[1];
    Type q = p[2];

    ctx.report("sigma", sigma);
    ctx.report("survival_probability", p_survive);
    ctx.report("q", q);

    Type target_sigma = Type(1.0);
    Type target_p = Type(0.5);
    Type target_q = Type(2.0);

    Type nll = Type(0.0);
    nll += Type(0.5) * (sigma - target_sigma) * (sigma - target_sigma);
    nll += Type(0.5) * (p_survive - target_p) * (p_survive - target_p);
    nll += Type(0.5) * (q - target_q) * (q - target_q);

    return nll;
  }

private:
  quadra::ParameterSet parameters_m;
};

int main() {
  SimpleTransformedModel model;
  auto result = quadra::evaluate_gradient(model, model.initials());

  std::cout << "objective = " << result.objective_value_m << "\n";
  std::cout << "gradient:\n";

  const auto names = model.parameter_names();
  for (size_t i = 0; i < result.gradient_m.size(); ++i) {
    std::cout << "  d/d " << names[i] << " = " << result.gradient_m[i] << "\n";
  }

  std::cout << "reports:\n";
  for (const auto &r : result.reports_m) {
    std::cout << "  " << r.name_m << " = " << r.value_m << "\n";
  }

  return 0;
}
