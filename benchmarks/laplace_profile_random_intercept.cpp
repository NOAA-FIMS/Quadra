#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "../core/laplace/laplace_profile.hpp"
#include "../core/model/quadra_model.hpp"

DECLARE_ADGRAPH();

class RandomInterceptModel : public quadra::QuadraModel<RandomInterceptModel> {
public:
  explicit RandomInterceptModel(std::vector<double> y) : y_m(std::move(y)) {
    parameters_m.add("mu", 0.0, quadra::ParameterTransform::Identity, false);
    parameters_m.add("u", 0.0, quadra::ParameterTransform::Identity, true);
  }

  std::vector<std::string> parameter_names_impl() const {
    return parameters_m.names();
  }

  const quadra::ParameterSet &parameters() const { return parameters_m; }

  template <typename Type>
  Type evaluate_impl(const std::vector<Type> &p,
                     quadra::ModelReportContext &) const {
    Type mu = p[0];
    Type u = p[1];
    Type nll = Type(0.0);

    for (double yi : y_m) {
      Type r = Type(yi) - (mu + u);
      nll += Type(0.5) * r * r;
    }

    nll += Type(0.5) * u * u;
    return nll;
  }

private:
  std::vector<double> y_m;
  quadra::ParameterSet parameters_m;
};

std::vector<double> simulate_data(size_t n) {
  std::mt19937 rng(1234);
  std::normal_distribution<double> dist(5.25, 1.0);

  std::vector<double> y(n);
  for (size_t i = 0; i < n; ++i)
    y[i] = dist(rng);
  return y;
}

int main() {
  std::cout << "\nQuadra Laplace component profiling benchmark\n\n";

  std::cout << std::setw(8) << "n" << std::setw(16) << "Newton ms"
            << std::setw(16) << "Hessian ms" << std::setw(18) << "LaplaceObj ms"
            << std::setw(18) << "ExactGrad ms" << std::setw(18)
            << "ExactLBFGS ms" << std::setw(14) << "LBFGS iter" << std::setw(16)
            << "theta_hat" << std::setw(16) << "u_hat" << "\n";

  std::cout << std::string(140, '-') << "\n";

  for (size_t n : std::vector<size_t>{10, 100, 1000, 5000, 10000}) {
    RandomInterceptModel model(simulate_data(n));

    auto components = quadra::profile_laplace_components(model, {4.7}, {0.0},
                                                         model.parameters());

    auto optimizer =
        quadra::profile_exact_lbfgs(model, {4.0}, {0.0}, model.parameters());

    std::cout << std::setw(8) << n << std::setw(16) << std::fixed
              << std::setprecision(3) << components.random_newton_ms_m
              << std::setw(16) << components.random_hessian_ms_m
              << std::setw(18) << components.laplace_objective_ms_m
              << std::setw(18) << components.exact_gradient_ms_m
              << std::setw(18) << optimizer.exact_lbfgs_ms_m << std::setw(14)
              << optimizer.exact_lbfgs_iterations_m << std::setw(16)
              << optimizer.theta_hat_m << std::setw(16) << optimizer.u_hat_m
              << "\n";
  }

  return 0;
}
