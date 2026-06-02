#include "../../core/laplace/laplace_exact_objective_gradient.hpp"
#include "../../include/quadra/quadra.hpp"

#include <Eigen/SparseCholesky>

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

DECLARE_ADGRAPH()

struct ExactGradientStateSpaceModel {
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
  std::ofstream csv("benchmarks/exact_laplace_gradient/"
                    "state_space_exact_gradient_benchmark.csv");

  if (!csv) {
    std::cerr << "failed to open exact gradient benchmark CSV\n";
    return 1;
  }

  csv << "n_state,n_random,total_ms,objective_ms,"
      << "tape_setup_ms,reverse_pass_ms,gradient_extract_ms,"
      << "total_gradient_ms,laplace_objective,gradient_norm,"
      << "hessian_nnz,hessian_density,factor_nnz,fill_ratio,"
      << "converged,logdet_ok,success\n";

  for (int n_state : {25, 50, 100, 250}) {
    ExactGradientStateSpaceModel model;
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

    quadra::LaplaceExactObjectiveGradientOptions options;
    options.include_logdet_derivative_m = false;

    quadra::LaplaceExactObjectiveGradientResult result;

    const double total_ms = elapsed_ms([&]() {
      result = quadra::evaluate_laplace_exact_objective_gradient(
          model, theta, random_initial, parameters, options);
    });

    const auto &H = result.hessian_random_m;

    const double hessian_density =
        H.rows() > 0 ? static_cast<double>(H.nonZeros()) /
                           static_cast<double>(H.rows() * H.cols())
                     : std::numeric_limits<double>::quiet_NaN();

    int factor_nnz = 0;
    double fill_ratio = std::numeric_limits<double>::quiet_NaN();

    if (result.converged_m && result.logdet_ok_m) {
      Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
      ldlt.compute(H);

      if (ldlt.info() == Eigen::Success) {
        const Eigen::SparseMatrix<double> L = ldlt.matrixL();

        factor_nnz = L.nonZeros();

        fill_ratio = H.nonZeros() > 0
                         ? static_cast<double>(factor_nnz) /
                               static_cast<double>(H.nonZeros())
                         : std::numeric_limits<double>::quiet_NaN();
      }
    }

    const bool success = result.converged_m && result.logdet_ok_m &&
                         std::isfinite(result.gradient_norm_m);

    csv << n_state << "," << n_state << "," << total_ms << ","
        << result.objective_ms_m << "," << result.tape_setup_ms_m << ","
        << result.reverse_pass_ms_m << "," << result.gradient_extract_ms_m
        << "," << result.total_gradient_ms_m << ","
        << result.laplace_objective_m << "," << result.gradient_norm_m << ","
        << H.nonZeros() << "," << hessian_density << "," << factor_nnz << ","
        << fill_ratio << "," << (result.converged_m ? 1 : 0) << ","
        << (result.logdet_ok_m ? 1 : 0) << "," << (success ? 1 : 0) << "\n";

    std::cout << "n_state=" << n_state << " total_ms=" << total_ms
              << " objective_ms=" << result.objective_ms_m
              << " reverse_ms=" << result.reverse_pass_ms_m
              << " extract_ms=" << result.gradient_extract_ms_m
              << " gradient_norm=" << result.gradient_norm_m
              << " hessian_nnz=" << H.nonZeros() << " fill_ratio=" << fill_ratio
              << " success=" << success << "\n";

    if (!success) {
      std::cerr << "exact gradient benchmark failed: " << result.message_m
                << "\n";
      return 1;
    }
  }

  return 0;
}
