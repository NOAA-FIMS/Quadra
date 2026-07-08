#pragma once

#include "../../../../../core/optimizer.hpp"

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace pifsc_bigeye_tuna {

inline std::vector<std::string> level9_fixed_effect_names()
{
  return {
      "log_r0",
      "log_m",
      "log_fbar",
      "log_q_longline",
      "log_q_purse_seine",
      "logit_sel_a50_longline",
      "log_sel_slope_longline"};
}

template <class Objective>
inline double level9_direct_joint_value_at_fixed(
    Objective &objective,
    const quadra::OptResult &fit,
    const std::vector<double> &theta)
{
  std::vector<double> full;
  full.reserve(theta.size() + fit.u_hat.size());

  for (const double v : theta)
    full.push_back(v);
  for (const double v : fit.u_hat)
    full.push_back(v);

  return objective(full);
}

template <class Objective>
inline Eigen::MatrixXd level9_direct_joint_fixed_hessian(
    Objective &objective,
    const quadra::OptResult &fit,
    double relative_step = 1.0e-3,
    double absolute_min_step = 1.0e-4)
{
  const std::size_t n = fit.par.size();
  if (n == 0)
    throw std::runtime_error("Level 22 geometry received empty fixed-effect vector");

  Eigen::MatrixXd h = Eigen::MatrixXd::Zero(
      static_cast<Eigen::Index>(n), static_cast<Eigen::Index>(n));

  const std::vector<double> base = fit.par;

  auto step_for = [&](double x) {
    return std::max(absolute_min_step,
                    relative_step * std::max(1.0, std::abs(x)));
  };

  const double f0 = level9_direct_joint_value_at_fixed(objective, fit, base);

  for (std::size_t i = 0; i < n; ++i)
  {
    const double hi = step_for(base[i]);

    std::vector<double> plus = base;
    std::vector<double> minus = base;
    plus[i] += hi;
    minus[i] -= hi;

    const double fp = level9_direct_joint_value_at_fixed(objective, fit, plus);
    const double fm = level9_direct_joint_value_at_fixed(objective, fit, minus);

    h(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(i)) =
        (fp - 2.0 * f0 + fm) / (hi * hi);

    for (std::size_t j = i + 1; j < n; ++j)
    {
      const double hj = step_for(base[j]);

      std::vector<double> pp = base;
      std::vector<double> pm = base;
      std::vector<double> mp = base;
      std::vector<double> mm = base;

      pp[i] += hi; pp[j] += hj;
      pm[i] += hi; pm[j] -= hj;
      mp[i] -= hi; mp[j] += hj;
      mm[i] -= hi; mm[j] -= hj;

      const double fpp = level9_direct_joint_value_at_fixed(objective, fit, pp);
      const double fpm = level9_direct_joint_value_at_fixed(objective, fit, pm);
      const double fmp = level9_direct_joint_value_at_fixed(objective, fit, mp);
      const double fmm = level9_direct_joint_value_at_fixed(objective, fit, mm);

      const double hij = (fpp - fpm - fmp + fmm) / (4.0 * hi * hj);
      h(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = hij;
      h(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(i)) = hij;
    }
  }

  return 0.5 * (h + h.transpose());
}

template <class Objective>
inline void write_level9_fixed_effect_geometry_report(
    const std::string &text_path,
    const std::string &csv_path,
    Objective &objective,
    const quadra::OptResult &fit)
{
  const auto all_names = level9_fixed_effect_names();
  if (fit.par.size() > all_names.size())
    throw std::runtime_error("Level 22 geometry has more fixed effects than names");

  std::vector<std::string> names(all_names.begin(),
                                 all_names.begin() + fit.par.size());

  const Eigen::MatrixXd h =
      level9_direct_joint_fixed_hessian(objective, fit);

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(h);
  const bool eig_ok = eig.info() == Eigen::Success;
  const auto evals = eig.eigenvalues();
  const auto evecs = eig.eigenvectors();

  Eigen::MatrixXd cov;
  bool covariance_ok = false;
  if (eig_ok && evals.size() > 0 && evals.minCoeff() > 0.0)
  {
    cov = h.inverse();
    covariance_ok = true;
  }

  Eigen::MatrixXd corr = Eigen::MatrixXd::Zero(h.rows(), h.cols());
  if (covariance_ok)
  {
    for (Eigen::Index i = 0; i < cov.rows(); ++i)
    {
      for (Eigen::Index j = 0; j < cov.cols(); ++j)
      {
        const double den = std::sqrt(std::abs(cov(i, i) * cov(j, j)));
        corr(i, j) =
            den > 0.0 ? cov(i, j) / den
                      : std::numeric_limits<double>::quiet_NaN();
      }
    }
  }

  std::ofstream txt(text_path);
  if (!txt)
    throw std::runtime_error("Could not open Level 22 fixed-effect geometry text: " +
                             text_path);

  std::ofstream csv(csv_path);
  if (!csv)
    throw std::runtime_error("Could not open Level 22 fixed-effect geometry CSV: " +
                             csv_path);

  txt << std::setprecision(15);
  csv << std::setprecision(15);

  txt << "Level 22 Fixed Effect Geometry Report\n";
  txt << "====================================\n\n";
  txt << "geometry_type:          direct joint objective at fixed u_hat\n";
  txt << "available:             yes\n";
  txt << "eigen_success:         " << (eig_ok ? "yes" : "no") << "\n";
  txt << "positive_definite:     "
      << (eig_ok && evals.minCoeff() > 0.0 ? "yes" : "no") << "\n";
  txt << "covariance_available:  " << (covariance_ok ? "yes" : "no") << "\n";
  if (eig_ok)
  {
    txt << "min_eigenvalue:        " << evals.minCoeff() << "\n";
    txt << "max_eigenvalue:        " << evals.maxCoeff() << "\n";
    txt << "condition_number_abs:  "
        << std::abs(evals.maxCoeff() / evals.minCoeff()) << "\n";
  }

  txt << "\nInterpretation note\n";
  txt << "-------------------\n";
  txt << "This is direct joint-objective geometry holding u_hat fixed, not fully\n";
  txt << "profiled-Laplace geometry. It is still useful for detecting local ridges\n";
  txt << "and parameter confounding around the returned Level 22 solution.\n\n";

  txt << "Fixed effects\n";
  txt << "-------------\n";
  txt << "index,name,value\n";
  for (std::size_t i = 0; i < fit.par.size(); ++i)
    txt << i << "," << names[i] << "," << fit.par[i] << "\n";

  txt << "\nEigenvalues\n";
  txt << "-----------\n";
  txt << "index,eigenvalue\n";
  if (eig_ok)
    for (Eigen::Index i = 0; i < evals.size(); ++i)
      txt << i << "," << evals[i] << "\n";

  txt << "\nWeak directions\n";
  txt << "---------------\n";
  if (eig_ok)
  {
    const Eigen::Index n_dirs = std::min<Eigen::Index>(3, evecs.cols());
    for (Eigen::Index d = 0; d < n_dirs; ++d)
    {
      std::vector<std::pair<double, std::size_t>> loadings;
      for (std::size_t i = 0; i < names.size(); ++i)
        loadings.push_back(
            {std::abs(evecs(static_cast<Eigen::Index>(i), d)), i});

      std::sort(loadings.begin(), loadings.end(),
                [](const auto &a, const auto &b) {
                  return a.first > b.first;
                });

      txt << "direction " << (d + 1) << ", eigenvalue = " << evals[d]
          << "\n";
      txt << "parameter,loading\n";
      for (const auto &entry : loadings)
      {
        const std::size_t i = entry.second;
        txt << names[i] << ","
            << evecs(static_cast<Eigen::Index>(i), d) << "\n";
      }
      txt << "\n";
    }
  }

  txt << "Fixed-effect correlation matrix\n";
  txt << "-------------------------------\n";
  txt << "row,col,value\n";
  if (covariance_ok)
  {
    for (std::size_t i = 0; i < names.size(); ++i)
      for (std::size_t j = 0; j < names.size(); ++j)
        txt << names[i] << "," << names[j] << ","
            << corr(static_cast<Eigen::Index>(i),
                    static_cast<Eigen::Index>(j))
            << "\n";
  }

  csv << "section,metric,target,value,extra\n";
  csv << "summary,geometry_type,,direct_joint_fixed_u_hat,\n";
  csv << "summary,eigen_success,," << (eig_ok ? "yes" : "no") << ",\n";
  csv << "summary,positive_definite,,"
      << (eig_ok && evals.minCoeff() > 0.0 ? "yes" : "no") << ",\n";
  csv << "summary,covariance_available,,"
      << (covariance_ok ? "yes" : "no") << ",\n";
  if (eig_ok)
  {
    csv << "summary,min_eigenvalue,," << evals.minCoeff() << ",\n";
    csv << "summary,max_eigenvalue,," << evals.maxCoeff() << ",\n";
    csv << "summary,condition_number_abs,,"
        << std::abs(evals.maxCoeff() / evals.minCoeff()) << ",\n";
  }

  for (std::size_t i = 0; i < fit.par.size(); ++i)
    csv << "fixed_effect,value," << names[i] << "," << fit.par[i]
        << ",index=" << i << "\n";

  if (eig_ok)
  {
    for (Eigen::Index i = 0; i < evals.size(); ++i)
      csv << "eigenvalue,value," << i << "," << evals[i] << ",\n";

    const Eigen::Index n_dirs = std::min<Eigen::Index>(3, evecs.cols());
    for (Eigen::Index d = 0; d < n_dirs; ++d)
      for (std::size_t i = 0; i < names.size(); ++i)
        csv << "weak_direction,loading,direction_" << (d + 1) << ","
            << evecs(static_cast<Eigen::Index>(i), d)
            << ",parameter=" << names[i] << ";eigenvalue=" << evals[d]
            << "\n";
  }

  if (covariance_ok)
    for (std::size_t i = 0; i < names.size(); ++i)
      for (std::size_t j = 0; j < names.size(); ++j)
        csv << "correlation,value," << names[i] << "_vs_" << names[j]
            << ","
            << corr(static_cast<Eigen::Index>(i),
                    static_cast<Eigen::Index>(j))
            << ",\n";
}

} // namespace pifsc_bigeye_tuna

using pifsc_bigeye_tuna::write_level9_fixed_effect_geometry_report;
