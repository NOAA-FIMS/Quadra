#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "../core/laplace/laplace_objective_cached.hpp"
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
  std::cout << "Testing cached Laplace objective\n";

  const int G = 10;
  const int m = 5;

  std::vector<double> y;
  std::vector<int> group;

  for (int g = 0; g < G; ++g) {
    for (int i = 0; i < m; ++i) {
      y.push_back(5.0);
      group.push_back(g);
    }
  }

  CorrelatedRandomInterceptModel model(y, group, G);

  std::vector<double> theta = {5.0};
  std::vector<double> u0(G, 0.0);

  auto uncached =
      quadra::evaluate_laplace_objective(model, theta, u0, model.parameters());

  quadra::CachedLaplaceObjectiveState cache_state;

  auto cached1 = quadra::evaluate_laplace_objective_cached(
      model, theta, u0, model.parameters(), cache_state);

  theta[0] = 5.1;

  auto cached2 = quadra::evaluate_laplace_objective_cached(
      model, theta, cached1.u_hat_m, model.parameters(), cache_state);

  std::cout << "uncached objective = " << uncached.laplace_objective_m << "\n";
  std::cout << "cached1 objective = " << cached1.laplace_objective_m << "\n";
  std::cout << "cached2 objective = " << cached2.laplace_objective_m << "\n";
  std::cout << "cached analyzed = " << cache_state.analyzed_m << "\n";
  std::cout << "cached nnz = " << cache_state.factorization_m.nonzeros()
            << "\n";

  if (!uncached.converged_m || !uncached.logdet_ok_m) {
    std::cerr << "FAIL: uncached objective failed\n";
    return 1;
  }

  if (!cached1.converged_m || !cached1.logdet_ok_m) {
    std::cerr << "FAIL: cached1 objective failed\n";
    return 1;
  }

  if (!cached2.converged_m || !cached2.logdet_ok_m) {
    std::cerr << "FAIL: cached2 objective failed\n";
    return 1;
  }

  if (!cache_state.analyzed_m) {
    std::cerr << "FAIL: cache was not analyzed\n";
    return 1;
  }

  if (std::abs(uncached.laplace_objective_m - cached1.laplace_objective_m) >
      1e-8) {
    std::cerr << "FAIL: cached and uncached objective mismatch\n";
    return 1;
  }

  std::cout << "PASS\n";
  return 0;
}
