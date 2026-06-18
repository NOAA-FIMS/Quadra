#pragma once

#include "../objective/red_snapper_quadra_objective.hpp"

#include "../../../../core/laplace/functional_analysis_report.hpp"
#include "../../../../core/laplace/laplace_structure_report.hpp"
#include "../../../../core/optimizer.hpp"

#include <Eigen/Dense>

#include <cmath>
#include <string>
#include <vector>

namespace sefsc_red_snapper {

inline std::size_t
red_snapper_max_fixed_gradient_index(const quadra::OptResult &fit) {
  std::size_t max_i = 0;
  double max_abs = -1.0;

  for (std::size_t i = 0; i < fit.fixed_gradient.size(); ++i) {
    const double a = std::abs(fit.fixed_gradient[i]);
    if (std::isfinite(a) && a > max_abs) {
      max_abs = a;
      max_i = i;
    }
  }

  return max_i;
}

template <class Objective>
Eigen::MatrixXd red_snapper_fd_huu(Objective &objective,
                                   const quadra::ParameterVector & /*params*/,
                                   const quadra::OptResult &fit,
                                   double rel_step = 1.0e-4) {
  const std::size_t n_fixed = fit.par.size();
  const std::size_t n_random = fit.u_hat.size();

  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(
      static_cast<Eigen::Index>(n_random), static_cast<Eigen::Index>(n_random));

  if (n_fixed == 0 || n_random == 0) {
    return H;
  }

  auto make_x = [&](const std::vector<double> &u) {
    std::vector<double> x;
    x.reserve(n_fixed + n_random);
    for (double v : fit.par) {
      x.push_back(v);
    }
    for (double v : u) {
      x.push_back(v);
    }
    return x;
  };

  auto eval_u = [&](const std::vector<double> &u) {
    return static_cast<double>(objective(make_x(u)));
  };

  const std::vector<double> u0 = fit.u_hat;
  const double f0 = eval_u(u0);

  std::vector<double> hi(n_random);
  for (std::size_t i = 0; i < n_random; ++i) {
    hi[i] = std::max(1.0e-5, rel_step * (1.0 + std::abs(u0[i])));
  }

  for (std::size_t i = 0; i < n_random; ++i) {
    auto up = u0;
    auto um = u0;
    up[i] += hi[i];
    um[i] -= hi[i];

    const double fp = eval_u(up);
    const double fm = eval_u(um);

    H(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(i)) =
        (fp - 2.0 * f0 + fm) / (hi[i] * hi[i]);
  }

  for (std::size_t i = 0; i < n_random; ++i) {
    for (std::size_t j = i + 1; j < n_random; ++j) {
      auto xpp = u0;
      auto xpm = u0;
      auto xmp = u0;
      auto xmm = u0;

      xpp[i] += hi[i];
      xpp[j] += hi[j];
      xpm[i] += hi[i];
      xpm[j] -= hi[j];
      xmp[i] -= hi[i];
      xmp[j] += hi[j];
      xmm[i] -= hi[i];
      xmm[j] -= hi[j];

      const double fpp = eval_u(xpp);
      const double fpm = eval_u(xpm);
      const double fmp = eval_u(xmp);
      const double fmm = eval_u(xmm);

      const double hij = (fpp - fpm - fmp + fmm) / (4.0 * hi[i] * hi[j]);

      H(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = hij;
      H(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(i)) = hij;
    }
  }

  return H;
}

template <class Objective>
void write_red_snapper_laplace_structure_report(
    const std::string &text_path, const std::string &csv_path,
    Objective &objective, const quadra::ParameterVector &params,
    const quadra::OptResult &fit, double nonzero_tol = 1.0e-8) {
  const Eigen::MatrixXd H = red_snapper_fd_huu(objective, params, fit);
  const auto report =
      quadra::summarize_laplace_hessian_structure(H, nonzero_tol);

  quadra::write_laplace_structure_report_text(report, text_path);
  quadra::write_laplace_structure_report_csv(report, csv_path);
}

template <class Objective>
void write_red_snapper_functional_analysis_report(
    const std::string &text_path, const std::string &csv_path,
    Objective &objective, const quadra::ParameterVector &params,
    const quadra::OptResult &fit, double nonzero_tol = 1.0e-8) {
  const Eigen::MatrixXd H = red_snapper_fd_huu(objective, params, fit);

  quadra::FunctionalOptimizationSummary opt;
  opt.objective_value = fit.value;
  opt.gradient_norm = fit.grad_norm;
  opt.iterations = fit.iterations;
  opt.converged = fit.converged;
  opt.message = fit.message;

  if (!fit.fixed_gradient.empty()) {
    const std::size_t max_i = red_snapper_max_fixed_gradient_index(fit);
    opt.max_gradient_parameter = (max_i < fit.fixed_gradient_names.size())
                                     ? fit.fixed_gradient_names[max_i]
                                     : ("fixed_" + std::to_string(max_i));
    opt.max_gradient_value = fit.fixed_gradient[max_i];
    opt.max_abs_gradient = std::abs(fit.fixed_gradient[max_i]);
  }

  std::vector<std::string> random_names;
  random_names.reserve(fit.u_hat.size());
  for (std::size_t i = 0; i < fit.u_hat.size(); ++i) {
    random_names.push_back("log_rec_dev_" + std::to_string(i + 1));
  }

  auto report = quadra::make_functional_analysis_report(
      opt, H, fit.u_hat, nonzero_tol, random_names);

  quadra::write_functional_analysis_report_text(report, text_path);
  quadra::write_functional_analysis_report_csv(report, csv_path);
}

} // namespace sefsc_red_snapper

using sefsc_red_snapper::red_snapper_fd_huu;
using sefsc_red_snapper::write_red_snapper_functional_analysis_report;
using sefsc_red_snapper::write_red_snapper_laplace_structure_report;
