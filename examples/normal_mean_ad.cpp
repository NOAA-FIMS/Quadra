#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "../core/autodiff/model_gradient.hpp"
#include "../core/model/quadra_model.hpp"

DECLARE_ADGRAPH();

class NormalMeanModel : public quadra::QuadraModel<NormalMeanModel> {
public:
  explicit NormalMeanModel(std::vector<double> y) : y_m(std::move(y)) {}

  std::vector<std::string> parameter_names_impl() const { return {"mu"}; }

  template <typename Type>
  Type evaluate_impl(const std::vector<Type> &parameters,
                     quadra::ModelReportContext &ctx) const {
    Type mu = parameters[0];
    Type nll = Type(0.0);

    for (double yi : y_m) {
      Type residual = Type(yi) - mu;
      nll += Type(0.5) * residual * residual;
    }

    ctx.report("mu", mu);
    ctx.adreport("mean_prediction", mu);
    return nll;
  }

private:
  std::vector<double> y_m;
};

int main() {
  NormalMeanModel model({4.8, 5.1, 5.0, 4.9, 5.2});
  std::vector<double> parameters = {4.7};

  auto result = quadra::evaluate_gradient(model, parameters);

  std::cout << "objective = " << result.objective_value_m << "\n";
  std::cout << "gradient[mu] = " << result.gradient_m[0] << "\n";

  for (const auto &report : result.reports_m) {
    std::cout << report.name_m << " = " << report.value_m;
    if (report.requires_se_m) {
      std::cout << " [ADREPORT]";
    }
    std::cout << "\n";
  }

  return 0;
}
