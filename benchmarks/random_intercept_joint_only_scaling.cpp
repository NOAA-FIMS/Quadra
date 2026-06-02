#include <chrono>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

#include "../core/laplace/joint_only_exact_lbfgs_optimizer.hpp"
#include "../core/laplace/laplace_exact_objective_lbfgs_optimizer.hpp"
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

template <typename F> double time_ms(F &&f) {
  auto start = std::chrono::high_resolution_clock::now();
  f();
  auto end = std::chrono::high_resolution_clock::now();
  return std::chrono::duration<double, std::milli>(end - start).count();
}

int main() {
  std::cout << "\nQuadra full Laplace exact vs joint-only exact benchmark\n\n";

  std::cout << std::setw(8) << "n" << std::setw(20) << "Full exact ms"
            << std::setw(14) << "Full iter" << std::setw(20) << "Joint only ms"
            << std::setw(14) << "Joint iter" << std::setw(16) << "theta_hat"
            << std::setw(16) << "u_hat" << "\n";

  std::cout << std::string(108, '-') << "\n";

  for (size_t n : std::vector<size_t>{10, 100, 1000, 5000, 10000}) {
    RandomInterceptModel model(simulate_data(n));

    quadra::LaplaceExactObjectiveLBFGSResult full_result;
    quadra::JointOnlyExactLBFGSResult joint_result;

    double full_ms = time_ms([&]() {
      full_result =
          quadra::optimize_laplace_fixed_effects_exact_objective_lbfgs(
              model, {4.0}, {0.0}, model.parameters());
    });

    double joint_ms = time_ms([&]() {
      joint_result = quadra::optimize_joint_only_exact_lbfgs(
          model, {4.0}, {0.0}, model.parameters());
    });

    std::cout << std::setw(8) << n << std::setw(20) << std::fixed
              << std::setprecision(3) << full_ms << std::setw(14)
              << full_result.iterations_m << std::setw(20) << joint_ms
              << std::setw(14) << joint_result.iterations_m << std::setw(16)
              << joint_result.theta_hat_m[0] << std::setw(16)
              << joint_result.u_hat_m[0] << "\n";
  }

  return 0;
}
