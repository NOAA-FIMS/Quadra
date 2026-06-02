#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "../core/laplace/laplace_evaluator.hpp"
#include "../core/model/quadra_model.hpp"

DECLARE_ADGRAPH();

class CorrelatedRandomInterceptModel
    : public quadra::QuadraModel<CorrelatedRandomInterceptModel> {
public:
  CorrelatedRandomInterceptModel(std::vector<double> y, std::vector<int> group,
                                 int n_groups)
      : y_m(std::move(y)), group_m(std::move(group)), n_groups_m(n_groups) {
    parameters_m.add("mu", 0.0, quadra::ParameterTransform::Identity, false);

    for (int g = 0; g < n_groups_m; ++g) {
      parameters_m.add("u_" + std::to_string(g), 0.0,
                       quadra::ParameterTransform::Identity, true);
    }
  }

  std::vector<std::string> parameter_names_impl() const {
    return parameters_m.names();
  }

  const quadra::ParameterSet &parameters() const { return parameters_m; }

  template <typename Type>
  Type evaluate_impl(const std::vector<Type> &p,
                     quadra::ModelReportContext &) const {
    Type mu = p[0];
    Type nll = Type(0.0);

    for (size_t i = 0; i < y_m.size(); ++i) {
      const int g = group_m[i];
      Type u = p[1 + g];

      Type r = Type(y_m[i]) - (mu + u);
      nll += Type(0.5) * r * r;
    }

    Type u0 = p[1];
    nll += Type(0.5) * u0 * u0;

    const double rho = 0.8;

    for (int g = 1; g < n_groups_m; ++g) {
      Type ug = p[1 + g];
      Type up = p[1 + g - 1];

      Type diff = ug - Type(rho) * up;
      nll += Type(0.5) * diff * diff;
    }

    return nll;
  }

private:
  std::vector<double> y_m;
  std::vector<int> group_m;
  int n_groups_m;
  quadra::ParameterSet parameters_m;
};

int main() {
  std::cout << "Testing LaplaceEvaluator\n";

  const int G = 10;
  const int m = 5;

  std::vector<double> y;
  std::vector<int> group;

  for (int g = 0; g < G; ++g) {
    for (int i = 0; i < m; ++i) {
      y.push_back(5.0 + 0.1 * std::sin(static_cast<double>(g)));
      group.push_back(g);
    }
  }

  CorrelatedRandomInterceptModel model(y, group, G);

  quadra::LaplaceEvaluator<CorrelatedRandomInterceptModel> evaluator(
      model, model.parameters(), std::vector<double>(G, 0.0));

  auto r1 = evaluator.evaluate({4.9});
  auto r2 = evaluator.evaluate({5.0});

  std::cout << "r1 converged = " << r1.converged_m << "\n";
  std::cout << "r2 converged = " << r2.converged_m << "\n";
  std::cout << "evaluations = " << evaluator.evaluations() << "\n";
  std::cout << "cache analyzed = " << evaluator.cache_state().analyzed_m
            << "\n";
  std::cout << "nnz = " << r2.hessian_random_m.nonZeros() << "\n";
  std::cout << "random start size = " << evaluator.random_start().size()
            << "\n";

  if (!r1.converged_m || !r1.logdet_ok_m)
    return 1;
  if (!r2.converged_m || !r2.logdet_ok_m)
    return 1;
  if (evaluator.evaluations() != 2)
    return 1;
  if (!evaluator.cache_state().analyzed_m)
    return 1;
  if (evaluator.random_start().size() != static_cast<size_t>(G))
    return 1;

  std::cout << "PASS\n";
  return 0;
}
