#include <cmath>
#include <iostream>
#include <string>
#include <vector>

#include "../core/laplace/laplace_objective.hpp"
#include "../core/model/quadra_model.hpp"

DECLARE_ADGRAPH();

class CorrelatedRandomInterceptModel
    : public quadra::QuadraModel<CorrelatedRandomInterceptModel> {
public:
  CorrelatedRandomInterceptModel(std::vector<double> y, std::vector<int> group,
                                 int n_groups, double rho, double lambda0,
                                 double lambda_diff)
      : y_m(std::move(y)), group_m(std::move(group)), n_groups_m(n_groups),
        rho_m(rho), lambda0_m(lambda0), lambda_diff_m(lambda_diff) {
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
    nll += Type(0.5) * Type(lambda0_m) * u0 * u0;

    for (int g = 1; g < n_groups_m; ++g) {
      Type ug = p[1 + g];
      Type up = p[1 + g - 1];

      Type diff = ug - Type(rho_m) * up;
      nll += Type(0.5) * Type(lambda_diff_m) * diff * diff;
    }

    return nll;
  }

private:
  std::vector<double> y_m;
  std::vector<int> group_m;
  int n_groups_m;
  double rho_m;
  double lambda0_m;
  double lambda_diff_m;
  quadra::ParameterSet parameters_m;
};

int main() {
  std::cout << "Testing correlated random-intercept Hessian structure\n";

  const int G = 5;
  const int m = 3;

  std::vector<double> y;
  std::vector<int> group;

  for (int g = 0; g < G; ++g) {
    for (int i = 0; i < m; ++i) {
      y.push_back(5.0);
      group.push_back(g);
    }
  }

  CorrelatedRandomInterceptModel model(y, group, G,
                                       0.8, // rho
                                       1.0, // lambda0
                                       1.0  // lambda_diff
  );

  std::vector<double> theta = {5.0};
  std::vector<double> u0(G, 0.0);

  auto result =
      quadra::evaluate_laplace_objective(model, theta, u0, model.parameters());

  std::cout << "converged = " << result.converged_m << "\n";
  std::cout << "nnz Huu = " << result.hessian_random_m.nonZeros() << "\n";
  std::cout << "logdet = " << result.log_det_hessian_m << "\n";

  // Tridiagonal symmetric matrix has:
  // diagonal: G
  // off-diagonal: 2*(G-1)
  // total: 3G - 2
  const int expected_nnz = 3 * G - 2;

  if (!result.converged_m) {
    std::cerr << "FAIL: Laplace objective did not converge\n";
    return 1;
  }

  if (result.hessian_random_m.nonZeros() != expected_nnz) {
    std::cerr << "FAIL: Huu nnz mismatch, expected " << expected_nnz << "\n";
    return 1;
  }

  std::cout << "PASS\n";
  return 0;
}
