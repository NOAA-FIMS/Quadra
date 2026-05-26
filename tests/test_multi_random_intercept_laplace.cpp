#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "../core/laplace/laplace_exact_objective_lbfgs_optimizer.hpp"
#include "../core/model/quadra_model.hpp"

DECLARE_ADGRAPH();

class MultiRandomInterceptModel
    : public quadra::QuadraModel<MultiRandomInterceptModel> {
public:
  MultiRandomInterceptModel(std::vector<double> y, std::vector<int> group,
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

    for (int g = 0; g < n_groups_m; ++g) {
      Type u = p[1 + g];
      nll += Type(0.5) * u * u;
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
  std::cout << "Testing multi-random-intercept Laplace stack\n";

  const int G = 3;
  const int m = 4;

  std::vector<double> y;
  std::vector<int> group;

  const double mu_true = 5.0;
  const double u_true[G] = {-0.2, 0.0, 0.2};

  for (int g = 0; g < G; ++g) {
    for (int i = 0; i < m; ++i) {
      y.push_back(mu_true + u_true[g]);
      group.push_back(g);
    }
  }

  MultiRandomInterceptModel model(y, group, G);

  std::vector<double> theta0 = {4.5};
  std::vector<double> u0(G, 0.0);

  auto result = quadra::optimize_laplace_fixed_effects_exact_objective_lbfgs(
      model, theta0, u0, model.parameters());

  std::cout << "converged = " << result.converged_m << "\n";
  std::cout << "message = " << result.message_m << "\n";
  std::cout << "iterations = " << result.iterations_m << "\n";
  std::cout << "theta_hat = " << result.theta_hat_m[0] << "\n";
  std::cout << "gradient_norm = " << result.gradient_norm_m << "\n";
  std::cout << "u_hat:";
  for (double u : result.u_hat_m) {
    std::cout << " " << u;
  }
  std::cout << "\n";

  if (!result.converged_m) {
    std::cerr << "FAIL: optimizer did not converge\n";
    return 1;
  }

  if (result.u_hat_m.size() != static_cast<size_t>(G)) {
    std::cerr << "FAIL: wrong number of random effects\n";
    return 1;
  }

  if (result.gradient_norm_m > 1e-5) {
    std::cerr << "FAIL: final gradient too large\n";
    return 1;
  }

  std::cout << "PASS\n";
  return 0;
}
