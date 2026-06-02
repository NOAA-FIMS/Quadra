#include "../../../include/quadra/quadra.hpp"

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

DECLARE_ADGRAPH()

struct RandomInterceptCompareModel {
  int n_obs = 1000;

  template <typename Context> void initialize(Context &) {}

  template <typename T, typename Context>
  T evaluate(const std::vector<T> &p, Context &) const {
    const T mu = p[0];
    const T log_sigma = p[1];
    const T u = p[2];

    const T sigma = exp(log_sigma);

    T nll = T(0.5) * u * u;

    for (int i = 0; i < n_obs; ++i) {
      const double yi = 1.0 + 0.10 * std::sin(0.01 * static_cast<double>(i));

      const T resid = (T(yi) - (mu + u)) / sigma;

      nll += T(0.5) * resid * resid + log_sigma;
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
  std::ofstream csv("benchmarks/comparisons/tmb_random_intercept/"
                    "comparison_outputs/quadra_random_intercept_compare.csv");

  if (!csv) {
    std::cerr << "failed to open Quadra comparison CSV\n";
    return 1;
  }

  csv << "engine,n_obs,n_random,objective_eval_ms,objective_reps,workspace_ms,"
         "implicit_derivatives_ms,factorization_ms,total_wall_ms,hessian_nnz,"
         "success\n";

  for (int n_obs : {100, 1000, 5000, 10000}) {
    RandomInterceptCompareModel model;
    model.n_obs = n_obs;

    quadra::ParameterSet parameters;
    parameters.add("mu", 0.0, quadra::ParameterTransform::Identity, false);
    parameters.add("log_sigma", std::log(0.5),
                   quadra::ParameterTransform::Identity, false);
    parameters.add("u", 0.0, quadra::ParameterTransform::Identity, true);

    std::vector<double> theta{0.0, std::log(0.5)};
    std::vector<double> random_initial{0.0};

    const int objective_reps = 1000;

    std::vector<double> full{theta[0], theta[1], random_initial[0]};

    volatile double objective_sink = 0.0;

    const double objective_eval_total_ms = elapsed_ms([&]() {
      for (int rep = 0; rep < objective_reps; ++rep) {
        quadra::ModelReportContext ctx;
        objective_sink += model.evaluate<double>(full, ctx);
      }
    });

    const double objective_eval_ms =
        objective_eval_total_ms / static_cast<double>(objective_reps);

    quadra::LaplaceImplicitWorkspace workspace;

    const double wall_ms = elapsed_ms([&]() {
      workspace = quadra::build_laplace_implicit_workspace(
          model, theta, random_initial, parameters);
    });

    csv << "quadra," << n_obs << "," << 1 << "," << objective_eval_ms << ","
        << objective_reps << "," << workspace.total_ms_m << ","
        << workspace.implicit_derivatives_ms_m << ","
        << workspace.factorization_ms_m << "," << wall_ms << ","
        << workspace.H_uu_m.nonZeros() << "," << (workspace.success_m ? 1 : 0)
        << "\n";

    std::cout << "quadra n_obs=" << n_obs
              << " objective_eval_ms=" << objective_eval_ms
              << " wall_ms=" << wall_ms
              << " workspace_ms=" << workspace.total_ms_m
              << " success=" << workspace.success_m << "\n";

    if (!workspace.success_m) {
      std::cerr << "workspace failed: " << workspace.message_m << "\n";
      return 1;
    }
  }

  return 0;
}
