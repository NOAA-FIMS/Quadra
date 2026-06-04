#include "state_space_surplus_production.hpp"

#include "core/laplace/model_analysis_report.hpp"

#include <Eigen/Sparse>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace ss = quadra_examples::state_space_surplus_production;

namespace {

struct Transformed {
  double r;
  double K;
  double q;
  double sigma_process;
  double sigma_index;
  double B0;
  double inv_var_process;
  double inv_var_index;
};

struct LogPredDerivatives {
  double log_pred;
  double d1;
  double d2;
};

Transformed transform(const ss::Parameters &par) {
  Transformed tr;
  tr.r = std::exp(par.log_r);
  tr.K = std::exp(par.log_K);
  tr.q = std::exp(par.log_q);
  tr.sigma_process = std::exp(par.log_sigma_process);
  tr.sigma_index = std::exp(par.log_sigma_index);

  const double B0_frac = 1.0 / (1.0 + std::exp(-par.logit_B0_frac));
  tr.B0 = B0_frac * tr.K;

  tr.inv_var_process = 1.0 / (tr.sigma_process * tr.sigma_process);
  tr.inv_var_index = 1.0 / (tr.sigma_index * tr.sigma_index);

  return tr;
}

LogPredDerivatives log_pred_derivatives(const double log_B,
                                        const double catch_t,
                                        const Transformed &tr) {
  const double B = std::exp(log_B);
  const double production = tr.r * B * (1.0 - B / tr.K);
  const double pred = std::max(B + production - catch_t, 1e-9);

  const double dB_dx = B;
  const double dprod_dB = tr.r * (1.0 - 2.0 * B / tr.K);
  const double dA_dB = 1.0 + dprod_dB;
  const double dA_dx = dA_dB * dB_dx;

  const double d2prod_dB2 = -2.0 * tr.r / tr.K;
  const double d2A_dx2 = d2prod_dB2 * B * B + dA_dB * B;

  LogPredDerivatives out;
  out.log_pred = std::log(pred);
  out.d1 = dA_dx / pred;
  out.d2 = d2A_dx2 / pred - (dA_dx * dA_dx) / (pred * pred);
  return out;
}

Eigen::VectorXd deterministic_initial_x(const ss::Data &data,
                                        const ss::Parameters &par) {
  const int n = static_cast<int>(data.catch_observed.size());
  const Transformed tr = transform(par);
  Eigen::VectorXd x(n - 1);

  double B = tr.B0;
  for (int t = 0; t < n - 1; ++t) {
    const double production = tr.r * B * (1.0 - B / tr.K);
    B = std::max(B + production -
                     data.catch_observed[static_cast<std::size_t>(t)],
                 1e-9);
    x[t] = std::log(B);
  }
  return x;
}

Eigen::SparseMatrix<double> analytic_hessian_xx(const ss::Data &data,
                                                const ss::Parameters &par,
                                                const Eigen::VectorXd &x) {
  const int n_state = static_cast<int>(x.size());
  const Transformed tr = transform(par);

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(static_cast<std::size_t>(3 * n_state));

  for (int k = 0; k < n_state; ++k) {
    double diag = 0.0;
    diag += tr.inv_var_process;
    diag += tr.inv_var_index;

    if (k + 1 < n_state) {
      const auto lp_next = log_pred_derivatives(
          x[k], data.catch_observed[static_cast<std::size_t>(k + 1)], tr);
      const double e_next = x[k + 1] - lp_next.log_pred;

      diag +=
          (lp_next.d1 * lp_next.d1 - e_next * lp_next.d2) * tr.inv_var_process;

      const double off = -lp_next.d1 * tr.inv_var_process;
      triplets.emplace_back(k, k + 1, off);
      triplets.emplace_back(k + 1, k, off);
    }

    triplets.emplace_back(k, k, diag);
  }

  Eigen::SparseMatrix<double> H(n_state, n_state);
  H.setFromTriplets(triplets.begin(), triplets.end());
  H.makeCompressed();
  return H;
}

} // namespace

int main() {
  const ss::Data data = ss::make_demo_data();
  const ss::Parameters par = ss::make_demo_parameters();
  const Eigen::VectorXd x = deterministic_initial_x(data, par);
  const Eigen::SparseMatrix<double> H = analytic_hessian_xx(data, par, x);

  quadra::laplace::StructureDetectorOptions opts;
  opts.prefer_dense_for_small_matrices = false;
  opts.dense_size_cutoff = 0;
  opts.banded_width_cutoff = 64;

  const auto report = quadra::laplace::analyze_hessian_structure(H, opts);
  report.Print(std::cout);
  std::cout << "\n";

  return 0;
}
