#include "../../include/quadra/quadra.hpp"

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

DECLARE_ADGRAPH()

struct RandomInterceptBenchmarkModel {
  int n_obs = 10;

  template <typename Context> void initialize(Context &) {}

  template <typename T, typename Context>
  T evaluate(const std::vector<T> &p, Context &) const {
    const T mu = p[0];
    const T log_sigma = p[1];
    const T u = p[2];

    const T sigma = exp(log_sigma);
    T nll = T(0.5) * u * u;

    for (int i = 0; i < n_obs; ++i) {
      const double yi = 1.0 + 0.1 * std::sin(static_cast<double>(i));
      const T resid = T(yi) - (mu + u);
      nll += T(0.5) * (resid / sigma) * (resid / sigma) + log_sigma;
    }

    nll += T(0.5) * (mu / T(10.0)) * (mu / T(10.0));
    nll += T(0.5) * (log_sigma / T(2.0)) * (log_sigma / T(2.0));

    return nll;
  }
};

template <typename Function> double elapsed_ms(Function &&f) {
  const auto start = std::chrono::steady_clock::now();
  f();
  const auto end = std::chrono::steady_clock::now();

  return std::chrono::duration<double, std::milli>(end - start).count();
}

int main() {
  std::ofstream csv(
      "benchmarks/random_intercept/random_intercept_benchmark.csv");

  if (!csv) {
    std::cerr << "failed to open benchmark CSV\n";
    return 1;
  }

  csv << "n_obs,workspace_ms,implicit_derivatives_ms,factorization_ms,total_"
         "wall_ms,success\n";

  for (int n_obs : {10, 100, 1000, 5000}) {
    RandomInterceptBenchmarkModel model;
    model.n_obs = n_obs;

    quadra::ParameterSet parameters;
    parameters.add("mu", 0.0, quadra::ParameterTransform::Identity, false);
    parameters.add("log_sigma", std::log(0.5),
                   quadra::ParameterTransform::Identity, false);
    parameters.add("u", 0.0, quadra::ParameterTransform::Identity, true);

    std::vector<double> theta{0.0, std::log(0.5)};
    std::vector<double> random_initial{0.0};

    quadra::LaplaceImplicitWorkspace workspace;

    const double wall_ms = elapsed_ms([&]() {
      workspace = quadra::build_laplace_implicit_workspace(
          model, theta, random_initial, parameters);
    });

    csv << n_obs << "," << workspace.total_ms_m << ","
        << workspace.implicit_derivatives_ms_m << ","
        << workspace.factorization_ms_m << "," << wall_ms << ","
        << (workspace.success_m ? 1 : 0) << "\n";

    std::cout << "n_obs=" << n_obs << " workspace_ms=" << workspace.total_ms_m
              << " wall_ms=" << wall_ms << " success=" << workspace.success_m
              << "\n";

    if (!workspace.success_m) {
      std::cerr << "workspace failed: " << workspace.message_m << "\n";
      return 1;
    }
  }

  std::cout
      << "wrote benchmarks/random_intercept/random_intercept_benchmark.csv\n";

  return 0;
}
