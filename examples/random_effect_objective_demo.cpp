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

    for (double yi : y_m) {
      Type r = Type(yi) - (mu + u);
      nll += Type(0.5) * r * r;
    }

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
  RandomInterceptModel model({4.8, 5.1, 5.0, 4.9, 5.2});

  const auto partition = quadra::partition_parameters(model.parameters());

  std::vector<double> theta = {4.7};
  std::vector<double> u = {0.0};

  auto result = quadra::evaluate_random_effect_objective_gradient(model, theta,
                                                                  u, partition);

  std::cout << "random-effect objective demo\n";
  std::cout << "objective = " << result.objective_value_m << "\n";
  std::cout << "gradient wrt u:\n";

  for (double gi : result.gradient_random_m) {
    std::cout << "  " << gi << "\n";
  }

  std::cout << "reports:\n";
  for (const auto &report : result.reports_m) {
    std::cout << "  " << report.name_m << " = " << report.value_m;

    if (report.requires_se_m) {
      std::cout << " [ADREPORT]";
    }

    std::cout << "\n";
  }

  return 0;
}
