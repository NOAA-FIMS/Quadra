#pragma once

#include "../../../../../core/optimizer.hpp"

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

namespace pifsc_bigeye_tuna {

struct FixedEffectWeakDirection {
  int direction = 0;
  double eigenvalue = std::numeric_limits<double>::quiet_NaN();
  std::vector<std::pair<std::string, double>> loadings;
};

struct FixedEffectGeometryReport {
  bool available = false;
  bool positive_definite = false;
  bool covariance_available = false;

  double min_eigenvalue = std::numeric_limits<double>::quiet_NaN();
  double max_eigenvalue = std::numeric_limits<double>::quiet_NaN();
  double condition_number_abs = std::numeric_limits<double>::quiet_NaN();

  Eigen::MatrixXd hessian;
  Eigen::MatrixXd covariance;
  Eigen::MatrixXd correlation;

  std::vector<std::string> names;
  std::vector<double> gradient;
  std::vector<double> eigenvalues;
  std::vector<FixedEffectWeakDirection> weak_directions;
};

template <class Objective>
double bigeye_profiled_laplace_value_at_fixed(Objective &objective,
                                              quadra::ParameterVector params,
                                              const std::vector<double> &theta,
                                              quadra::LaplaceOptions opts) {
  for (std::size_t i = 0; i < theta.size() && i < params.params.size(); ++i) {
    params.params[i].value = theta[i];
  }

  std::vector<int> fixed_idx;
  std::vector<int> random_idx;

  for (std::size_t i = 0; i < params.params.size(); ++i) {
    if (params.params[i].is_random)
      random_idx.push_back(static_cast<int>(i));
    else
      fixed_idx.push_back(static_cast<int>(i));
  }

  Eigen::VectorXd x(static_cast<Eigen::Index>(fixed_idx.size()));
  for (std::size_t i = 0; i < fixed_idx.size(); ++i) {
    x[static_cast<Eigen::Index>(i)] = params.params[fixed_idx[i]].value;
  }

  had::ADGraph graph;
  const auto u_hat = quadra::solve_random_effects_laplace(
      objective, params, x, fixed_idx, random_idx, graph);

  const auto res = quadra::laplace_eval_at_u_star(
      objective, params, fixed_idx, random_idx, x, u_hat, graph, opts);

  return res.value;
}

template <class Objective>
Eigen::MatrixXd bigeye_fd_fixed_hessian(Objective &objective,
                                        const quadra::ParameterVector &params,
                                        const quadra::OptResult &fit,
                                        quadra::LaplaceOptions opts,
                                        double rel_step = 1.0e-4) {
  const std::size_t n = fit.par.size();
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(static_cast<Eigen::Index>(n),
                                            static_cast<Eigen::Index>(n));

  if (n == 0)
    return H;

  const std::vector<double> theta0 = fit.par;

  auto value_at = [&](const std::vector<double> &theta) {
    return bigeye_profiled_laplace_value_at_fixed(objective, params, theta,
                                                  opts);
  };

  const double f0 = value_at(theta0);

  std::vector<double> h(n);
  for (std::size_t i = 0; i < n; ++i) {
    h[i] = std::max(1.0e-5, rel_step * (1.0 + std::abs(theta0[i])));
  }

  for (std::size_t i = 0; i < n; ++i) {
    auto xp = theta0;
    auto xm = theta0;
    xp[i] += h[i];
    xm[i] -= h[i];

    const double fp = value_at(xp);
    const double fm = value_at(xm);

    H(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(i)) =
        (fp - 2.0 * f0 + fm) / (h[i] * h[i]);
  }

  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = i + 1; j < n; ++j) {
      auto xpp = theta0;
      auto xpm = theta0;
      auto xmp = theta0;
      auto xmm = theta0;

      xpp[i] += h[i];
      xpp[j] += h[j];
      xpm[i] += h[i];
      xpm[j] -= h[j];
      xmp[i] -= h[i];
      xmp[j] += h[j];
      xmm[i] -= h[i];
      xmm[j] -= h[j];

      const double fpp = value_at(xpp);
      const double fpm = value_at(xpm);
      const double fmp = value_at(xmp);
      const double fmm = value_at(xmm);

      const double hij = (fpp - fpm - fmp + fmm) / (4.0 * h[i] * h[j]);

      H(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) = hij;
      H(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(i)) = hij;
    }
  }

  return H;
}

inline std::vector<std::string>
bigeye_fixed_effect_names(const quadra::ParameterVector &params,
                          const quadra::OptResult &fit) {
  if (fit.fixed_gradient_names.size() == fit.par.size())
    return fit.fixed_gradient_names;

  std::vector<std::string> names;
  for (const auto &p : params.params) {
    if (!p.is_random)
      names.push_back(p.name);
  }

  while (names.size() < fit.par.size())
    names.push_back("fixed_" + std::to_string(names.size()));

  return names;
}

inline Eigen::MatrixXd covariance_to_correlation(const Eigen::MatrixXd &cov) {
  Eigen::MatrixXd cor = Eigen::MatrixXd::Zero(cov.rows(), cov.cols());

  for (Eigen::Index i = 0; i < cov.rows(); ++i) {
    for (Eigen::Index j = 0; j < cov.cols(); ++j) {
      const double vi = cov(i, i);
      const double vj = cov(j, j);
      const double denom = std::sqrt(std::max(0.0, vi) * std::max(0.0, vj));

      cor(i, j) = denom > 0.0 ? cov(i, j) / denom
                              : std::numeric_limits<double>::quiet_NaN();
    }
  }

  return cor;
}

template <class Objective>
FixedEffectGeometryReport make_bigeye_fixed_effect_geometry_report(
    Objective &objective, const quadra::ParameterVector &params,
    const quadra::OptResult &fit, quadra::LaplaceOptions opts,
    int n_weak_directions = 3) {
  FixedEffectGeometryReport report;
  report.names = bigeye_fixed_effect_names(params, fit);
  report.gradient = fit.fixed_gradient;

  if (fit.par.empty())
    return report;

  report.hessian = bigeye_fd_fixed_hessian(objective, params, fit, opts);
  report.available = true;

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(report.hessian);
  if (es.info() != Eigen::Success)
    return report;

  const auto evals = es.eigenvalues();
  const auto evecs = es.eigenvectors();

  for (Eigen::Index i = 0; i < evals.size(); ++i)
    report.eigenvalues.push_back(evals[i]);

  report.min_eigenvalue =
      evals.size() > 0 ? evals[0] : std::numeric_limits<double>::quiet_NaN();
  report.max_eigenvalue = evals.size() > 0
                              ? evals[evals.size() - 1]
                              : std::numeric_limits<double>::quiet_NaN();

  report.positive_definite = report.min_eigenvalue > 0.0;
  report.condition_number_abs =
      std::abs(report.max_eigenvalue) /
      std::max(std::abs(report.min_eigenvalue), 1.0e-300);

  if (report.positive_definite) {
    report.covariance = report.hessian.inverse();
    report.correlation = covariance_to_correlation(report.covariance);
    report.covariance_available = true;
  }

  const int nd =
      std::min<int>(n_weak_directions, static_cast<int>(report.names.size()));

  for (int d = 0; d < nd; ++d) {
    FixedEffectWeakDirection wd;
    wd.direction = d + 1;
    wd.eigenvalue = evals[d];

    for (Eigen::Index i = 0; i < evecs.rows(); ++i) {
      const std::size_t ui = static_cast<std::size_t>(i);
      const std::string name = ui < report.names.size()
                                   ? report.names[ui]
                                   : ("fixed_" + std::to_string(ui));
      wd.loadings.push_back({name, evecs(i, d)});
    }

    std::sort(wd.loadings.begin(), wd.loadings.end(),
              [](const auto &a, const auto &b) {
                return std::abs(a.second) > std::abs(b.second);
              });

    report.weak_directions.push_back(wd);
  }

  return report;
}

inline void
write_bigeye_fixed_effect_geometry_text(const FixedEffectGeometryReport &report,
                                        const std::string &path) {
  std::ofstream out(path);
  if (!out)
    throw std::runtime_error("Could not open fixed-effect geometry report: " +
                             path);

  out << "Fixed Effect Geometry Report\n";
  out << "============================\n\n";

  out << std::setprecision(15);
  out << "available:             " << (report.available ? "yes" : "no") << "\n";
  out << "positive_definite:     " << (report.positive_definite ? "yes" : "no")
      << "\n";
  out << "covariance_available:  "
      << (report.covariance_available ? "yes" : "no") << "\n";
  out << "min_eigenvalue:        " << report.min_eigenvalue << "\n";
  out << "max_eigenvalue:        " << report.max_eigenvalue << "\n";
  out << "condition_number_abs:  " << report.condition_number_abs << "\n\n";

  out << "Fixed effects\n";
  out << "-------------\n";
  out << "index,name,gradient\n";
  for (std::size_t i = 0; i < report.names.size(); ++i) {
    const double g = i < report.gradient.size()
                         ? report.gradient[i]
                         : std::numeric_limits<double>::quiet_NaN();
    out << i << "," << report.names[i] << "," << g << "\n";
  }

  out << "\nEigenvalues\n";
  out << "-----------\n";
  out << "index,eigenvalue\n";
  for (std::size_t i = 0; i < report.eigenvalues.size(); ++i)
    out << i << "," << report.eigenvalues[i] << "\n";

  out << "\nWeak directions\n";
  out << "---------------\n";
  for (const auto &wd : report.weak_directions) {
    out << "direction " << wd.direction << ", eigenvalue = " << wd.eigenvalue
        << "\n";
    out << "parameter,loading\n";
    for (const auto &x : wd.loadings)
      out << x.first << "," << x.second << "\n";
    out << "\n";
  }

  if (report.covariance_available) {
    out << "Fixed-effect correlation matrix\n";
    out << "-------------------------------\n";
    out << "row,col,value\n";
    for (Eigen::Index i = 0; i < report.correlation.rows(); ++i) {
      for (Eigen::Index j = 0; j < report.correlation.cols(); ++j) {
        out << report.names[static_cast<std::size_t>(i)] << ","
            << report.names[static_cast<std::size_t>(j)] << ","
            << report.correlation(i, j) << "\n";
      }
    }
  }

  out << "\nInterpretation notes\n";
  out << "--------------------\n";
  out << "Small eigenvalues identify weakly curved fixed-effect directions.\n";
  out << "Large loadings in the same weak direction indicate parameter "
         "combinations\n";
  out << "that may be difficult for the data to separate.\n";
}

inline void
write_bigeye_fixed_effect_geometry_csv(const FixedEffectGeometryReport &report,
                                       const std::string &path) {
  std::ofstream out(path);
  if (!out)
    throw std::runtime_error("Could not open fixed-effect geometry CSV: " +
                             path);

  out << "section,metric,target,value,extra\n";
  out << std::setprecision(15);

  out << "summary,available,," << (report.available ? "yes" : "no") << ",\n";
  out << "summary,positive_definite,,"
      << (report.positive_definite ? "yes" : "no") << ",\n";
  out << "summary,covariance_available,,"
      << (report.covariance_available ? "yes" : "no") << ",\n";
  out << "summary,min_eigenvalue,," << report.min_eigenvalue << ",\n";
  out << "summary,max_eigenvalue,," << report.max_eigenvalue << ",\n";
  out << "summary,condition_number_abs,," << report.condition_number_abs
      << ",\n";

  for (std::size_t i = 0; i < report.names.size(); ++i) {
    const double g = i < report.gradient.size()
                         ? report.gradient[i]
                         : std::numeric_limits<double>::quiet_NaN();
    out << "fixed_effect,gradient," << report.names[i] << "," << g
        << ",index=" << i << "\n";
  }

  for (std::size_t i = 0; i < report.eigenvalues.size(); ++i)
    out << "eigenvalue,value," << i << "," << report.eigenvalues[i] << ",\n";

  for (const auto &wd : report.weak_directions) {
    for (const auto &x : wd.loadings) {
      out << "weak_direction,loading,direction_" << wd.direction << ","
          << x.second << ",parameter=" << x.first
          << ";eigenvalue=" << wd.eigenvalue << "\n";
    }
  }

  if (report.covariance_available) {
    for (Eigen::Index i = 0; i < report.correlation.rows(); ++i) {
      for (Eigen::Index j = 0; j < report.correlation.cols(); ++j) {
        out << "correlation,value," << report.names[static_cast<std::size_t>(i)]
            << "_vs_" << report.names[static_cast<std::size_t>(j)] << ","
            << report.correlation(i, j) << ",\n";
      }
    }
  }
}

template <class Objective>
void write_bigeye_fixed_effect_geometry_report(
    const std::string &text_path, const std::string &csv_path,
    Objective &objective, const quadra::ParameterVector &params,
    const quadra::OptResult &fit, quadra::LaplaceOptions opts) {
  const auto report =
      make_bigeye_fixed_effect_geometry_report(objective, params, fit, opts);
  write_bigeye_fixed_effect_geometry_text(report, text_path);
  write_bigeye_fixed_effect_geometry_csv(report, csv_path);
}

} // namespace pifsc_bigeye_tuna

using pifsc_bigeye_tuna::FixedEffectGeometryReport;
using pifsc_bigeye_tuna::write_bigeye_fixed_effect_geometry_report;
