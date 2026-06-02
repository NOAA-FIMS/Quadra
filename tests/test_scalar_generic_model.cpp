#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "../core/autodiff/model_gradient.hpp"
#include "../core/model/quadra_model.hpp"
#include "../core/runtime/model_runner.hpp"

DECLARE_ADGRAPH();

class NormalMeanGenericModel
    : public quadra::QuadraModel<NormalMeanGenericModel> {
public:
  explicit NormalMeanGenericModel(std::vector<double> y) : y_m(std::move(y)) {}

  std::vector<std::string> parameter_names_impl() const { return {"mu"}; }

  template <typename Type>
  Type evaluate_impl(const std::vector<Type> &parameters,
                     quadra::ModelReportContext &ctx) const {
    Type mu = parameters[0];
    Type nll = Type(0.0);

    for (double yi : y_m) {
      Type r = Type(yi) - mu;
      nll += Type(0.5) * r * r;
    }

    ctx.report("mu", mu);
    ctx.report("n", static_cast<double>(y_m.size()));
    ctx.adreport("mean_prediction", mu);

    return nll;
  }

private:
  std::vector<double> y_m;
};

int main() {
  std::cout << "Testing scalar-generic model evaluation\n";

  NormalMeanGenericModel model({4.8, 5.1, 5.0, 4.9, 5.2});
  std::vector<double> parameters = {5.0};

  auto value_result = quadra::evaluate_once(model, parameters);
  auto grad_result = quadra::evaluate_gradient(model, parameters);

  const double expected_value = 0.05;
  const double expected_gradient = 0.0;

  std::cout << "objective = " << value_result.objective_value_m << "\n";
  std::cout << "gradient  = " << grad_result.gradient_m[0] << "\n";
  std::cout << "reports   = " << grad_result.reports_m.size() << "\n";

  if (std::abs(value_result.objective_value_m - expected_value) > 1e-12) {
    std::cerr << "FAIL: objective mismatch\n";
    return 1;
  }

  if (std::abs(grad_result.gradient_m[0] - expected_gradient) > 1e-12) {
    std::cerr << "FAIL: gradient mismatch\n";
    return 1;
  }

  if (grad_result.reports_m.size() != 3) {
    std::cerr << "FAIL: expected three reports\n";
    return 1;
  }

  std::cout << "PASS\n";
  return 0;
}
