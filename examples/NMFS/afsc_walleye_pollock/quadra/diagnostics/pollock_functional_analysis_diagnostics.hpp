#pragma once

#include "../diagnostics/pollock_fixed_effect_diagnostics.hpp"
#include "../diagnostics/pollock_fixed_hessian_diagnostics.hpp"
#include "../diagnostics/pollock_huu_diagnostics.hpp"
#include "../model/pollock_laplace_helpers.hpp"
#include "../model/pollock_model.hpp"

#include "../../../../../core/laplace/functional_analysis_report.hpp"
#include "../../../../../core/laplace/laplace_structure_report.hpp"
#include "../../../../../core/optimizer.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace pollock_example {

// Shared low-level model-evaluation/index helpers.
// These must appear before FD/Huu and higher-level diagnostics.

// Low-level dependencies used by the functional-analysis diagnostics.
// Keep these before the higher-level diagnostics that call them.

void pollock_write_huu_sparsity(const std::string &path, PollockModel &model,
                                quadra::ParameterVector params,
                                const quadra::OptResult &fit,
                                double tol = 1.0e-8) {
  std::ofstream out(path);
  out << "i,j,value,abs_value\n";
  out << std::setprecision(15);

  Eigen::MatrixXd dense = pollock_fd_huu(model, params, fit);

  for (Eigen::Index i = 0; i < dense.rows(); ++i) {
    for (Eigen::Index j = 0; j < dense.cols(); ++j) {
      const double v = dense(i, j);
      if (std::abs(v) > tol) {
        out << (i + 1) << "," << (j + 1) << "," << v << "," << std::abs(v)
            << "\n";
      }
    }
  }
}

void pollock_write_huu_band_summary(const std::string &path,
                                    PollockModel &model,
                                    quadra::ParameterVector params,
                                    const quadra::OptResult &fit,
                                    double tol = 1.0e-8) {
  Eigen::MatrixXd H = pollock_fd_huu(model, params, fit);

  std::ofstream out(path);
  out << "band_distance,count,nonzero_count,mean_abs,max_abs,sum_abs,share_sum_"
         "abs,cumulative_share_sum_abs\n";
  out << std::setprecision(15);

  if (H.rows() == 0) {
    return;
  }

  const Eigen::Index n = H.rows();
  std::vector<double> sum_abs(static_cast<std::size_t>(n), 0.0);
  std::vector<double> max_abs(static_cast<std::size_t>(n), 0.0);
  std::vector<std::size_t> count(static_cast<std::size_t>(n), 0);
  std::vector<std::size_t> nonzero_count(static_cast<std::size_t>(n), 0);

  double total_abs = 0.0;

  // Use upper triangle including diagonal so each symmetric pair is counted
  // once.
  for (Eigen::Index i = 0; i < n; ++i) {
    for (Eigen::Index j = i; j < n; ++j) {
      const std::size_t d = static_cast<std::size_t>(j - i);
      const double av = std::abs(H(i, j));

      count[d] += 1;
      sum_abs[d] += av;
      max_abs[d] = std::max(max_abs[d], av);
      total_abs += av;

      if (av > tol) {
        nonzero_count[d] += 1;
      }
    }
  }

  double cumulative = 0.0;
  for (std::size_t d = 0; d < static_cast<std::size_t>(n); ++d) {
    const double mean_abs =
        count[d] > 0 ? sum_abs[d] / static_cast<double>(count[d]) : 0.0;
    const double share = total_abs > 0.0 ? sum_abs[d] / total_abs : 0.0;
    cumulative += share;

    out << d << "," << count[d] << "," << nonzero_count[d] << "," << mean_abs
        << "," << max_abs[d] << "," << sum_abs[d] << "," << share << ","
        << cumulative << "\n";
  }
}

void pollock_write_huu_bandlimit_diagnostic(const std::string &path,
                                            PollockModel &model,
                                            quadra::ParameterVector params,
                                            const quadra::OptResult &fit) {
  Eigen::MatrixXd H = pollock_fd_huu(model, params, fit);

  std::ofstream out(path);
  out << "band_width,kept_entries,total_entries,kept_entry_share,"
         "retained_abs_share,relative_frobenius_error,"
         "min_eigenvalue,max_eigenvalue,positive_definite,condition_number_"
         "abs\n";
  out << std::setprecision(15);

  if (H.rows() == 0) {
    return;
  }

  const Eigen::Index n = H.rows();
  const double full_abs_sum = H.cwiseAbs().sum();
  const double full_frob = H.norm();

  const std::vector<int> bands = {0, 1, 2, 3, 5, 10, 20};

  for (const int bw_raw : bands) {
    const Eigen::Index bw =
        std::min<Eigen::Index>(static_cast<Eigen::Index>(bw_raw), n - 1);

    Eigen::MatrixXd Hb = Eigen::MatrixXd::Zero(n, n);
    std::size_t kept_entries = 0;

    for (Eigen::Index i = 0; i < n; ++i) {
      for (Eigen::Index j = 0; j < n; ++j) {
        if (std::abs(i - j) <= bw) {
          Hb(i, j) = H(i, j);
          ++kept_entries;
        }
      }
    }

    const double retained_abs_share =
        full_abs_sum > 0.0 ? Hb.cwiseAbs().sum() / full_abs_sum : 0.0;
    const double rel_frob_error =
        full_frob > 0.0 ? (H - Hb).norm() / full_frob : 0.0;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(Hb);
    const bool eig_ok = eig.info() == Eigen::Success;
    double min_eval = std::numeric_limits<double>::quiet_NaN();
    double max_eval = std::numeric_limits<double>::quiet_NaN();
    bool pd = false;
    double cond = std::numeric_limits<double>::quiet_NaN();

    if (eig_ok && eig.eigenvalues().size() > 0) {
      min_eval = eig.eigenvalues().minCoeff();
      max_eval = eig.eigenvalues().maxCoeff();
      pd = min_eval > 0.0;
      cond = std::abs(max_eval) / std::max(std::abs(min_eval), 1.0e-300);
    }

    out << bw << "," << kept_entries << "," << static_cast<std::size_t>(n * n)
        << "," << static_cast<double>(kept_entries) / static_cast<double>(n * n)
        << "," << retained_abs_share << "," << rel_frob_error << "," << min_eval
        << "," << max_eval << "," << (pd ? "yes" : "no") << "," << cond << "\n";
  }
}

void pollock_write_huu_threshold_diagnostic(const std::string &path,
                                            PollockModel &model,
                                            quadra::ParameterVector params,
                                            const quadra::OptResult &fit) {
  Eigen::MatrixXd H = pollock_fd_huu(model, params, fit);

  std::ofstream out(path);
  out << "threshold_type,threshold,absolute_threshold,"
         "kept_entries,total_entries,kept_entry_share,"
         "retained_abs_share,relative_frobenius_error,"
         "min_eigenvalue,max_eigenvalue,positive_definite,condition_number_"
         "abs\n";
  out << std::setprecision(15);

  if (H.rows() == 0) {
    return;
  }

  const Eigen::Index n = H.rows();
  const double full_abs_sum = H.cwiseAbs().sum();
  const double full_frob = H.norm();
  const double max_abs = H.cwiseAbs().maxCoeff();

  struct ThresholdSpec {
    const char *type;
    double threshold;
    double absolute_threshold;
  };

  std::vector<ThresholdSpec> specs;

  for (double t : {1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5, 1.0e-6}) {
    specs.push_back({"absolute", t, t});
  }

  for (double r : {1.0e-1, 1.0e-2, 1.0e-3, 1.0e-4, 1.0e-5}) {
    specs.push_back({"relative_to_max_abs", r, r * max_abs});
  }

  for (const auto &spec : specs) {
    Eigen::MatrixXd Ht = Eigen::MatrixXd::Zero(n, n);
    std::size_t kept_entries = 0;

    for (Eigen::Index i = 0; i < n; ++i) {
      for (Eigen::Index j = 0; j < n; ++j) {
        const double v = H(i, j);
        if (std::abs(v) >= spec.absolute_threshold) {
          Ht(i, j) = v;
          ++kept_entries;
        }
      }
    }

    const double retained_abs_share =
        full_abs_sum > 0.0 ? Ht.cwiseAbs().sum() / full_abs_sum : 0.0;
    const double rel_frob_error =
        full_frob > 0.0 ? (H - Ht).norm() / full_frob : 0.0;

    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(Ht);
    const bool eig_ok = eig.info() == Eigen::Success;
    double min_eval = std::numeric_limits<double>::quiet_NaN();
    double max_eval = std::numeric_limits<double>::quiet_NaN();
    bool pd = false;
    double cond = std::numeric_limits<double>::quiet_NaN();

    if (eig_ok && eig.eigenvalues().size() > 0) {
      min_eval = eig.eigenvalues().minCoeff();
      max_eval = eig.eigenvalues().maxCoeff();
      pd = min_eval > 0.0;
      cond = std::abs(max_eval) / std::max(std::abs(min_eval), 1.0e-300);
    }

    out << spec.type << "," << spec.threshold << "," << spec.absolute_threshold
        << "," << kept_entries << "," << static_cast<std::size_t>(n * n) << ","
        << static_cast<double>(kept_entries) / static_cast<double>(n * n) << ","
        << retained_abs_share << "," << rel_frob_error << "," << min_eval << ","
        << max_eval << "," << (pd ? "yes" : "no") << "," << cond << "\n";
  }
}

void pollock_write_laplace_structure_report(const std::string &path,
                                            PollockModel &model,
                                            quadra::ParameterVector params,
                                            const quadra::OptResult &fit,
                                            double nonzero_tol = 1.0e-8) {
  const Eigen::MatrixXd H = pollock_fd_huu(model, params, fit);
  const auto report =
      quadra::summarize_laplace_hessian_structure(H, nonzero_tol);

  quadra::write_laplace_structure_report_text(report, path);
  quadra::write_laplace_structure_report_csv(
      report, "examples/NMFS/afsc_walleye_pollock/outputs/"
              "walleye_pollock_laplace_structure_report.csv");
}

quadra::FunctionalGradientVolatilitySummary
pollock_compute_gradient_volatility_fd(PollockModel &model,
                                       quadra::ParameterVector params,
                                       const quadra::OptResult &fit,
                                       double perturbation_scale = 1.0e-5,
                                       double fd_step = 1.0e-5) {
  quadra::FunctionalGradientVolatilitySummary empty;
  if (fit.fixed_gradient.empty()) {
    return empty;
  }

  Eigen::VectorXd x0 = fit.x;
  if (x0.size() == 0) {
    x0 = pollock_fixed_values(params);
  }

  if (x0.size() == 0 ||
      x0.size() != static_cast<Eigen::Index>(fit.fixed_gradient.size())) {
    return empty;
  }

  quadra::LaplaceOptions opts = quadra::default_laplace_options();

  std::vector<std::vector<double>> gradient_samples;
  gradient_samples.reserve(static_cast<std::size_t>(x0.size() * 2));

  for (Eigen::Index j = 0; j < x0.size(); ++j) {
    const double dx = perturbation_scale * std::max(1.0, std::abs(x0(j)));

    for (double sign : {-1.0, 1.0}) {
      Eigen::VectorXd xp = x0;
      xp(j) += sign * dx;
      gradient_samples.push_back(pollock_profile_gradient_fd_at_x(
          model, params, xp, fit.u_hat, opts, fd_step));
    }
  }

  return quadra::summarize_gradient_volatility(
      gradient_samples, fit.fixed_gradient, fit.fixed_gradient_names,
      perturbation_scale);
}

std::vector<int> pollock_parameter_geometry_fixed_indices(
    const quadra::ParameterVector &params) {
  std::vector<int> out;
  for (std::size_t i = 0; i < params.params.size(); ++i) {
    if (!params.params[i].is_random) {
      out.push_back(static_cast<int>(i));
    }
  }
  return out;
}

std::vector<int> pollock_parameter_geometry_random_indices(
    const quadra::ParameterVector &params) {
  std::vector<int> out;
  for (std::size_t i = 0; i < params.params.size(); ++i) {
    if (params.params[i].is_random) {
      out.push_back(static_cast<int>(i));
    }
  }
  return out;
}

Eigen::MatrixXd pollock_parameter_geometry_fd_fixed_hessian(
    PollockModel &model, quadra::ParameterVector params,
    const quadra::OptResult &fit, const quadra::LaplaceOptions &opts,
    double fd_step = 1.0e-4) {
  Eigen::VectorXd x0 = fit.x;
  if (x0.size() == 0) {
    // Backward-compatible fallback. Prefer fit.x when available.
    const auto fixed_idx = pollock_parameter_geometry_fixed_indices(params);
    x0.resize(static_cast<Eigen::Index>(fixed_idx.size()));
    for (std::size_t i = 0; i < fixed_idx.size(); ++i) {
      x0(static_cast<Eigen::Index>(i)) =
          params.params[static_cast<std::size_t>(fixed_idx[i])].value;
    }
  }

  const Eigen::Index n = x0.size();
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(n, n);

  if (n == 0) {
    return H;
  }

  const auto fixed_idx = pollock_parameter_geometry_fixed_indices(params);
  const auto random_idx = pollock_parameter_geometry_random_indices(params);

  auto eval = [&](const Eigen::VectorXd &x_eval) -> double {
    had::ADGraph graph;
    auto res = quadra::laplace_eval_at_u_star(
        model, params, fixed_idx, random_idx, x_eval, fit.u_hat, graph, opts);
    return res.value;
  };

  for (Eigen::Index i = 0; i < n; ++i) {
    const double hi = fd_step * std::max(1.0, std::abs(x0(i)));

    // Diagonal second derivative.
    {
      Eigen::VectorXd xp = x0;
      Eigen::VectorXd xm = x0;
      xp(i) += hi;
      xm(i) -= hi;

      const double f0 = eval(x0);
      const double fp = eval(xp);
      const double fm = eval(xm);
      H(i, i) = (fp - 2.0 * f0 + fm) / (hi * hi);
    }

    // Mixed second derivatives.
    for (Eigen::Index j = i + 1; j < n; ++j) {
      const double hj = fd_step * std::max(1.0, std::abs(x0(j)));

      Eigen::VectorXd xpp = x0;
      Eigen::VectorXd xpm = x0;
      Eigen::VectorXd xmp = x0;
      Eigen::VectorXd xmm = x0;

      xpp(i) += hi;
      xpp(j) += hj;
      xpm(i) += hi;
      xpm(j) -= hj;
      xmp(i) -= hi;
      xmp(j) += hj;
      xmm(i) -= hi;
      xmm(j) -= hj;

      const double fpp = eval(xpp);
      const double fpm = eval(xpm);
      const double fmp = eval(xmp);
      const double fmm = eval(xmm);

      const double hij = (fpp - fpm - fmp + fmm) / (4.0 * hi * hj);
      H(i, j) = hij;
      H(j, i) = hij;
    }
  }

  return H;
}

void pollock_write_functional_analysis_report(const std::string &text_path,
                                              const std::string &csv_path,
                                              PollockModel &model,
                                              quadra::ParameterVector params,
                                              const quadra::OptResult &fit,
                                              double nonzero_tol = 1.0e-8) {
  const Eigen::MatrixXd H = pollock_fd_huu(model, params, fit);

  quadra::FunctionalOptimizationSummary opt;
  opt.objective_value = fit.value;
  opt.gradient_norm = fit.grad_norm;
  opt.iterations = fit.iterations;
  opt.converged = fit.converged;
  opt.message = fit.message;

  if (!fit.fixed_gradient.empty()) {
    const std::size_t max_i = max_fixed_gradient_index(fit);
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

#ifdef WALLEYE_POLLOCK_PARAMETER_GEOMETRY
  {
    quadra::LaplaceOptions hess_opts = quadra::default_laplace_options();
    const Eigen::MatrixXd Hxx = pollock_parameter_geometry_fd_fixed_hessian(
        model, params, fit, hess_opts);

    report.parameter_geometry = quadra::summarize_parameter_geometry(
        Hxx, fit.fixed_gradient, fit.fixed_gradient_names);
  }
#endif

#ifdef WALLEYE_POLLOCK_GRADIENT_VOLATILITY
  {
    report.gradient_volatility = pollock_compute_gradient_volatility_fd(
        model, params, fit, 1.0e-5, 1.0e-5);
  }
#endif

  quadra::write_functional_analysis_report_text(report, text_path);
  quadra::write_functional_analysis_report_csv(report, csv_path);
}

} // namespace pollock_example

// Compatibility aliases for current walleye_pollock.cpp call sites.
using pollock_example::pollock_compute_gradient_volatility_fd;
using pollock_example::pollock_parameter_geometry_fd_fixed_hessian;
using pollock_example::pollock_parameter_geometry_fixed_indices;
using pollock_example::pollock_parameter_geometry_random_indices;
using pollock_example::pollock_write_functional_analysis_report;
using pollock_example::pollock_write_huu_band_summary;
using pollock_example::pollock_write_huu_bandlimit_diagnostic;
using pollock_example::pollock_write_huu_sparsity;
using pollock_example::pollock_write_huu_threshold_diagnostic;
using pollock_example::pollock_write_laplace_structure_report;
