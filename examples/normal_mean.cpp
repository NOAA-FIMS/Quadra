#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

#include "../core/model/quadra_model.hpp"
#include "../core/runtime/model_runner.hpp"

class NormalMeanModel : public quadra::QuadraModelBase {
public:
  explicit NormalMeanModel(std::vector<double> y) : y_m(std::move(y)) {}

  std::vector<std::string> parameter_names() const override { return {"mu"}; }

  double evaluate(const std::vector<double> &parameters,
                  quadra::ModelContext &ctx) override {
    double mu = parameters[0];
    double nll = 0.0;

    for (double yi : y_m) {
      double r = yi - mu;
      nll += 0.5 * r * r;
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
  NormalMeanModel model({4.8, 5.1, 5.0, 4.9, 5.2});

  quadra::ModelRunner runner(model);

  std::vector<double> parameters = {5.0};

  quadra::ModelResult result = runner.evaluate_once(parameters);

  std::cout << "Objective value: " << result.objective_value_m << "\\n";

  std::cout << "Reports:\\n";
  for (const auto &r : result.reports_m) {
    std::cout << "  " << r.name_m << " = " << r.value_m;

    if (r.requires_se_m) {
      std::cout << "  [ADREPORT]";
    }

    std::cout << "\\n";
  }

  return 0;
}
