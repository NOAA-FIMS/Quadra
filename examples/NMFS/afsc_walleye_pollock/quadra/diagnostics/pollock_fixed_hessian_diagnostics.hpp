#pragma once

#include "../model/pollock_model.hpp"
#include "../model/pollock_laplace_helpers.hpp"

#include "../../../../../core/optimizer.hpp"
#include "../../../../../core/laplace/laplace_structure_report.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <exception>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>

#ifdef WALLEYE_POLLOCK_FIXED_HESSIAN_DIAGNOSTICS

namespace pollock_example {

double pollock_profile_objective_at_fixed(
  PollockModel &model,
  quadra::ParameterVector params,
  const std::vector<int> &fixed_idx,
  const std::vector<int> &random_idx,
  const Eigen::VectorXd &x,
  const quadra::LaplaceOptions &opts) {
for (std::size_t k = 0; k < fixed_idx.size(); ++k) {
  params.params[static_cast<std::size_t>(fixed_idx[k])].value =
      x[static_cast<Eigen::Index>(k)];
}

had::ADGraph graph;

if (random_idx.empty()) {
  std::vector<double> p_double;
  p_double.reserve(static_cast<std::size_t>(params.size()));
  for (int i = 0; i < params.size(); ++i) {
    p_double.emplace_back(
        params.params[static_cast<std::size_t>(i)].value);
  }
  return model(p_double);
}

const auto u_hat = quadra::solve_random_effects_laplace(
    model, params, x, fixed_idx, random_idx, graph);
const auto res = quadra::laplace_eval_at_u_star(
    model, params, fixed_idx, random_idx, x, u_hat, graph, opts);
return res.value;
}

void write_fixed_hessian_diagnostics(
  const std::string &summary_path,
  const std::string &matrix_path,
  PollockModel &model,
  const quadra::ParameterVector &params_in,
  const quadra::OptResult &fit,
  const quadra::LaplaceOptions &opts) {
std::ofstream summary(summary_path);
summary << std::setprecision(15);
summary << "field,value\n";

quadra::ParameterVector params = params_in;
const auto fixed_idx = quadra::build_fixed_index(params);
const auto random_idx = quadra::build_random_index(params);

const Eigen::Index n = static_cast<Eigen::Index>(fixed_idx.size());
summary << "fixed_effects," << n << "\n";

if (n == 0 || fit.par.size() != static_cast<std::size_t>(n)) {
  summary << "available,no\n";
  summary << "reason,missing fixed-effect vector\n";
  return;
}

try {
  Eigen::VectorXd x(n);
  for (Eigen::Index i = 0; i < n; ++i) {
    x[i] = fit.par[static_cast<std::size_t>(i)];
  }

  const double eps = 1.0e-4;
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(n, n);

  const double f0 = pollock_profile_objective_at_fixed(
      model, params, fixed_idx, random_idx, x, opts);

  for (Eigen::Index i = 0; i < n; ++i) {
    Eigen::VectorXd xp = x;
    Eigen::VectorXd xm = x;
    xp[i] += eps;
    xm[i] -= eps;

    const double fp = pollock_profile_objective_at_fixed(
        model, params, fixed_idx, random_idx, xp, opts);
    const double fm = pollock_profile_objective_at_fixed(
        model, params, fixed_idx, random_idx, xm, opts);

    H(i, i) = (fp - 2.0 * f0 + fm) / (eps * eps);

    for (Eigen::Index j = i + 1; j < n; ++j) {
      Eigen::VectorXd xpp = x;
      Eigen::VectorXd xpm = x;
      Eigen::VectorXd xmp = x;
      Eigen::VectorXd xmm = x;

      xpp[i] += eps;
      xpp[j] += eps;
      xpm[i] += eps;
      xpm[j] -= eps;
      xmp[i] -= eps;
      xmp[j] += eps;
      xmm[i] -= eps;
      xmm[j] -= eps;

      const double fpp = pollock_profile_objective_at_fixed(
          model, params, fixed_idx, random_idx, xpp, opts);
      const double fpm = pollock_profile_objective_at_fixed(
          model, params, fixed_idx, random_idx, xpm, opts);
      const double fmp = pollock_profile_objective_at_fixed(
          model, params, fixed_idx, random_idx, xmp, opts);
      const double fmm = pollock_profile_objective_at_fixed(
          model, params, fixed_idx, random_idx, xmm, opts);

      const double hij = (fpp - fpm - fmp + fmm) / (4.0 * eps * eps);
      H(i, j) = hij;
      H(j, i) = hij;
    }
  }

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(H);

  summary << "available,yes\n";
  summary << "fd_step," << eps << "\n";
  summary << "profile_objective," << f0 << "\n";
  summary << "min_diagonal," << H.diagonal().minCoeff() << "\n";
  summary << "max_diagonal," << H.diagonal().maxCoeff() << "\n";

  if (es.info() == Eigen::Success) {
    const auto evals = es.eigenvalues();
    summary << "eigen_success,yes\n";
    summary << "min_eigenvalue," << evals.minCoeff() << "\n";
    summary << "max_eigenvalue," << evals.maxCoeff() << "\n";
    summary << "positive_definite,"
            << (evals.minCoeff() > 0.0 ? "yes" : "no") << "\n";

    if (std::abs(evals.minCoeff()) > 0.0) {
      summary << "condition_number_abs,"
              << std::abs(evals.maxCoeff()) / std::abs(evals.minCoeff())
              << "\n";
    } else {
      summary << "condition_number_abs,inf\n";
    }

    summary << "eigenvalues";
    for (Eigen::Index i = 0; i < evals.size(); ++i) {
      summary << (i == 0 ? "," : ";") << evals[i];
    }
    summary << "\n";
  } else {
    summary << "eigen_success,no\n";
  }

  std::ofstream mat(matrix_path);
  mat << "parameter";
  for (std::size_t j = 0; j < fit.fixed_gradient_names.size(); ++j) {
    mat << "," << fit.fixed_gradient_names[j];
  }
  mat << "\n";

  for (Eigen::Index i = 0; i < n; ++i) {
    const std::string row_name =
        (static_cast<std::size_t>(i) < fit.fixed_gradient_names.size())
            ? fit.fixed_gradient_names[static_cast<std::size_t>(i)]
            : ("fixed_" + std::to_string(i));
    mat << row_name;
    for (Eigen::Index j = 0; j < n; ++j) {
      mat << "," << std::setprecision(15) << H(i, j);
    }
    mat << "\n";
  }
} catch (const std::exception &e) {
  summary << "available,no\n";
  summary << "reason," << e.what() << "\n";
}
}

}  // namespace pollock_example

using pollock_example::pollock_profile_objective_at_fixed;
using pollock_example::write_fixed_hessian_diagnostics;

#endif  // WALLEYE_POLLOCK_FIXED_HESSIAN_DIAGNOSTICS
