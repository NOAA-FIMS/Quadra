#include "../../core/laplace/laplace_implicit_workspace.hpp"
#include "../../core/laplace/sparse_factorization_cache.hpp"
#include "../../include/quadra/quadra.hpp"

#include <Eigen/SparseCholesky>

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

DECLARE_ADGRAPH()

struct FactorReuseStateSpaceModel {
  int n_state = 100;

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
      "benchmarks/exact_laplace_gradient/factorization_reuse_benchmark.csv");

  if (!csv) {
    std::cerr << "failed to open factorization reuse benchmark CSV\n";
    return 1;
  }

  csv << "n_state,iteration,fresh_ms,reuse_ms,reuse_ratio,"
      << "hessian_nnz,factor_nnz,fill_ratio,success\n";

  for (int n_state : {25, 50, 100, 250, 500}) {
    FactorReuseStateSpaceModel model;
    model.n_state = n_state;

    quadra::ParameterSet parameters;
    parameters.add("mu", 1.0, quadra::ParameterTransform::Identity, false);
    parameters.add("log_sigma_obs", std::log(0.4),
                   quadra::ParameterTransform::Identity, false);
    parameters.add("log_sigma_rw", std::log(0.3),
                   quadra::ParameterTransform::Identity, false);

    std::vector<double> theta{1.0, std::log(0.4), std::log(0.3)};
    std::vector<double> random_initial;

    for (int t = 0; t < n_state; ++t) {
      parameters.add("x_" + std::to_string(t), 0.0,
                     quadra::ParameterTransform::Identity, true);

      random_initial.push_back(0.0);
    }

    auto workspace = quadra::build_laplace_implicit_workspace(
        model, theta, random_initial, parameters);

    if (!workspace.success_m) {
      std::cerr << "workspace failed: " << workspace.message_m << "\n";
      return 1;
    }

    const auto &H = workspace.H_uu_m;

    int factor_nnz = 0;
    double fill_ratio = 0.0;

    {
      Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
      ldlt.compute(H);

      if (ldlt.info() != Eigen::Success) {
        std::cerr << "initial LDLT failed\n";
        return 1;
      }

      const Eigen::SparseMatrix<double> L = ldlt.matrixL();
      factor_nnz = L.nonZeros();
      fill_ratio =
          static_cast<double>(factor_nnz) / static_cast<double>(H.nonZeros());
    }

    quadra::SparseLDLTFactorizationCache cache;
    cache.analyze_pattern(H);

    for (int iter = 0; iter < 20; ++iter) {
      double fresh_ms = 0.0;
      double reuse_ms = 0.0;

      bool fresh_ok = true;
      bool reuse_ok = true;

      fresh_ms = elapsed_ms([&]() {
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> ldlt;
        ldlt.compute(H);

        if (ldlt.info() != Eigen::Success) {
          fresh_ok = false;
        }
      });

      reuse_ms = elapsed_ms([&]() {
        try {
          cache.factorize(H);
        } catch (...) {
          reuse_ok = false;
        }
      });

      const double reuse_ratio = fresh_ms > 0.0 ? reuse_ms / fresh_ms : 0.0;

      const bool success = fresh_ok && reuse_ok;

      csv << n_state << "," << iter << "," << fresh_ms << "," << reuse_ms << ","
          << reuse_ratio << "," << H.nonZeros() << "," << factor_nnz << ","
          << fill_ratio << "," << (success ? 1 : 0) << "\n";

      if (!success) {
        std::cerr << "factorization reuse benchmark failed\n";
        return 1;
      }
    }

    std::cout << "n_state=" << n_state << " hessian_nnz=" << H.nonZeros()
              << " factor_nnz=" << factor_nnz << " fill_ratio=" << fill_ratio
              << "\n";
  }

  return 0;
}
