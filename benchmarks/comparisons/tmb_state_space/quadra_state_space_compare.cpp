#include "../../../include/quadra/quadra.hpp"

#include <Eigen/SparseCholesky>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

DECLARE_ADGRAPH()

struct GaussianRandomWalkCompareModel {
  int n_state = 25;

  template <typename Context> void initialize(Context &) {}

  template <typename T, typename Context>
  T evaluate(const std::vector<T> &p, Context &) const {
    const T mu = p[0];
    const T log_sigma_obs = p[1];
    const T log_sigma_rw = p[2];

    const T sigma_obs = exp(log_sigma_obs);
    const T sigma_rw = exp(log_sigma_rw);

    T nll = T(0.0);

    for (int t = 0; t < n_state; ++t) {
      const T x_t = p[3 + t];

      const double y_t = 1.0 + 0.10 * std::sin(0.03 * static_cast<double>(t)) +
                         0.05 * std::cos(0.11 * static_cast<double>(t));

      const T obs_resid = (T(y_t) - (mu + x_t)) / sigma_obs;

      nll += T(0.5) * obs_resid * obs_resid + log_sigma_obs;

      if (t == 0) {
        const T z = x_t / sigma_rw;
        nll += T(0.5) * z * z + log_sigma_rw;
      } else {
        const T x_prev = p[3 + t - 1];
        const T z = (x_t - x_prev) / sigma_rw;
        nll += T(0.5) * z * z + log_sigma_rw;
      }
    }

    nll += T(0.001) * mu * mu;
    nll += T(0.001) * log_sigma_obs * log_sigma_obs;
    nll += T(0.001) * log_sigma_rw * log_sigma_rw;

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
  std::ofstream csv("benchmarks/comparisons/tmb_state_space/comparison_outputs/"
                    "quadra_state_space_compare.csv");

  if (!csv) {
    std::cerr << "failed to open Quadra state-space comparison CSV\n";
    return 1;
  }

  csv << "engine,model,n_state,n_random,objective_eval_ms,objective_reps,"
         "workspace_ms,implicit_derivatives_ms,factorization_ms,total_wall_ms,"
         "hessian_nnz,hessian_density,factor_nnz,fill_ratio,success\n";

  for (int n_state : {25, 50, 100, 250}) {
    GaussianRandomWalkCompareModel model;
    model.n_state = n_state;

    quadra::ParameterSet parameters;
    parameters.add("mu", 1.0, quadra::ParameterTransform::Identity, false);
    parameters.add("log_sigma_obs", std::log(0.4),
                   quadra::ParameterTransform::Identity, false);
    parameters.add("log_sigma_rw", std::log(0.3),
                   quadra::ParameterTransform::Identity, false);

    std::vector<double> theta{1.0, std::log(0.4), std::log(0.3)};
    std::vector<double> random_initial;
    random_initial.reserve(static_cast<std::size_t>(n_state));

    for (int t = 0; t < n_state; ++t) {
      parameters.add("x_" + std::to_string(t), 0.0,
                     quadra::ParameterTransform::Identity, true);

      random_initial.push_back(0.0);
    }

    std::vector<double> full;
    full.reserve(static_cast<std::size_t>(3 + n_state));
    full.push_back(theta[0]);
    full.push_back(theta[1]);
    full.push_back(theta[2]);

    for (int t = 0; t < n_state; ++t) {
      full.push_back(0.0);
    }

    const int objective_reps = n_state <= 50 ? 1000 : 200;

    volatile double objective_sink = 0.0;

    const double objective_total_ms = elapsed_ms([&]() {
      for (int rep = 0; rep < objective_reps; ++rep) {
        quadra::ModelReportContext ctx;
        objective_sink += model.evaluate<double>(full, ctx);
      }
    });

    const double objective_eval_ms =
        objective_total_ms / static_cast<double>(objective_reps);

    quadra::LaplaceImplicitWorkspace workspace;

    const double wall_ms = elapsed_ms([&]() {
      workspace = quadra::build_laplace_implicit_workspace(
          model, theta, random_initial, parameters);
    });

    const double hessian_density =
        workspace.H_uu_m.rows() > 0
            ? static_cast<double>(workspace.H_uu_m.nonZeros()) /
                  static_cast<double>(workspace.H_uu_m.rows() *
                                      workspace.H_uu_m.cols())
            : std::numeric_limits<double>::quiet_NaN();

    int factor_nnz = 0;
    double fill_ratio = std::numeric_limits<double>::quiet_NaN();

    if (workspace.success_m) {
      Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
      ldlt.compute(workspace.H_uu_m);

      if (ldlt.info() == Eigen::Success) {
        const Eigen::SparseMatrix<double> L = ldlt.matrixL();

        factor_nnz = L.nonZeros();

        fill_ratio = workspace.H_uu_m.nonZeros() > 0
                         ? static_cast<double>(factor_nnz) /
                               static_cast<double>(workspace.H_uu_m.nonZeros())
                         : std::numeric_limits<double>::quiet_NaN();
      }
    }

    csv << "quadra,state_space," << n_state << "," << n_state << ","
        << objective_eval_ms << "," << objective_reps << ","
        << workspace.total_ms_m << "," << workspace.implicit_derivatives_ms_m
        << "," << workspace.factorization_ms_m << "," << wall_ms << ","
        << workspace.H_uu_m.nonZeros() << "," << hessian_density << ","
        << factor_nnz << "," << fill_ratio << ","
        << (workspace.success_m ? 1 : 0) << "\n";

    std::cout << "quadra state_space n_state=" << n_state
              << " objective_eval_ms=" << objective_eval_ms
              << " workspace_ms=" << workspace.total_ms_m
              << " nnz=" << workspace.H_uu_m.nonZeros()
              << " density=" << hessian_density << " factor_nnz=" << factor_nnz
              << " fill_ratio=" << fill_ratio
              << " success=" << workspace.success_m << "\n";

    if (!workspace.success_m) {
      std::cerr << "workspace failed: " << workspace.message_m << "\n";
      return 1;
    }
  }

  return 0;
}
