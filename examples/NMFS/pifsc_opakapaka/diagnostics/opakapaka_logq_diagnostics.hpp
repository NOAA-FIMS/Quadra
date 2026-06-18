#pragma once

#include "../quadra/opakapaka_model.hpp"

#include "../../../../core/uncertainty/reporting.hpp"
#include "../../../../core/uncertainty/selected_inverse_diagonal.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <string>
#include <vector>

namespace opakapaka_example {

struct LogQUncertaintyReport {
  double objective = std::numeric_limits<double>::quiet_NaN();
  double fd_step = std::numeric_limits<double>::quiet_NaN();
  double fd_gradient = std::numeric_limits<double>::quiet_NaN();
  double fd_hessian = std::numeric_limits<double>::quiet_NaN();
  double covariance_log_q = std::numeric_limits<double>::quiet_NaN();
  double se_log_q = std::numeric_limits<double>::quiet_NaN();
  double log_q = std::numeric_limits<double>::quiet_NaN();
  double q = std::numeric_limits<double>::quiet_NaN();
  double se_q = std::numeric_limits<double>::quiet_NaN();
  double log_q_lwr_95 = std::numeric_limits<double>::quiet_NaN();
  double log_q_upr_95 = std::numeric_limits<double>::quiet_NaN();
  double q_lwr_95 = std::numeric_limits<double>::quiet_NaN();
  double q_upr_95 = std::numeric_limits<double>::quiet_NaN();
};

template <class Model>
LogQUncertaintyReport
compute_log_q_uncertainty_report(Model &model, quadra::ParameterVector &params,
                                 quadra::LaplaceOptions &opts,
                                 const quadra::OptResult &fit) {
  LogQUncertaintyReport out;
  if (fit.par.size() != 1)
    return out;

  const std::vector<int> fixed_idx = {0};
  std::vector<int> random_idx;
  for (std::size_t i = 1; i < params.size(); ++i) {
    random_idx.push_back(static_cast<int>(i));
  }

  auto eval_at = [&](double theta) {
    auto tmp = params;
    tmp.params.at(0).value = theta;
    Eigen::VectorXd x(1);
    x[0] = theta;
    had::ADGraph graph;
    auto u_hat = quadra::solve_random_effects_laplace(model, tmp, x, fixed_idx,
                                                      random_idx, graph);
    auto res = quadra::laplace_eval_at_u_star(model, tmp, fixed_idx, random_idx,
                                              x, u_hat, graph, opts);
    return res.value;
  };

  out.objective = fit.value;
  out.log_q = fit.par.at(0);
  out.q = std::exp(out.log_q);
  out.fd_step = std::max(1.0e-5, 1.0e-4 * (1.0 + std::abs(out.log_q)));

  const double fm = eval_at(out.log_q - out.fd_step);
  const double fp = eval_at(out.log_q + out.fd_step);
  if (!std::isfinite(fm) || !std::isfinite(fp) || !std::isfinite(out.objective))
    return out;

  out.fd_gradient = (fp - fm) / (2.0 * out.fd_step);
  out.fd_hessian =
      (fp - 2.0 * out.objective + fm) / (out.fd_step * out.fd_step);

  if (std::isfinite(out.fd_hessian) && out.fd_hessian > 0.0) {
    out.covariance_log_q = 1.0 / out.fd_hessian;
    out.se_log_q = std::sqrt(out.covariance_log_q);
    out.se_q = out.q * out.se_log_q;
    out.log_q_lwr_95 = out.log_q - 1.96 * out.se_log_q;
    out.log_q_upr_95 = out.log_q + 1.96 * out.se_log_q;
    out.q_lwr_95 = std::exp(out.log_q_lwr_95);
    out.q_upr_95 = std::exp(out.log_q_upr_95);
  }
  return out;
}

inline void write_uncertainty_summary_csv(const std::string &path,
                                          const LogQUncertaintyReport &u) {
  std::ofstream out(path);
  out << "field,value\n";
  out << "objective," << u.objective << "\n";
  out << "fd_step," << u.fd_step << "\n";
  out << "fd_gradient_log_q," << u.fd_gradient << "\n";
  out << "fd_hessian_log_q," << u.fd_hessian << "\n";
  out << "covariance_log_q," << u.covariance_log_q << "\n";
  out << "se_log_q," << u.se_log_q << "\n";
  out << "se_q," << u.se_q << "\n";
  out << "hessian_positive," << (u.fd_hessian > 0.0 ? "yes" : "no") << "\n";
}

inline void write_covariance_matrix_csv(const std::string &path,
                                        const LogQUncertaintyReport &u) {
  std::ofstream out(path);
  out << "row,col,value\n";
  out << "log_q,log_q," << u.covariance_log_q << "\n";
}

inline void write_correlation_matrix_csv(const std::string &path) {
  std::ofstream out(path);
  out << "row,col,value\n";
  out << "log_q,log_q,1\n";
}

inline void write_standard_errors_csv(const std::string &path,
                                      const LogQUncertaintyReport &u) {
  std::ofstream out(path);
  out << "parameter,scale,estimate,se\n";
  out << "log_q,log," << u.log_q << "," << u.se_log_q << "\n";
  out << "q,natural," << u.q << "," << u.se_q << "\n";
}

inline void write_confidence_intervals_csv(const std::string &path,
                                           const LogQUncertaintyReport &u) {
  std::ofstream out(path);
  out << "parameter,scale,estimate,se,lwr_95,upr_95\n";
  out << "log_q,log," << u.log_q << "," << u.se_log_q << "," << u.log_q_lwr_95
      << "," << u.log_q_upr_95 << "\n";
  out << "q,natural," << u.q << "," << u.se_q << "," << u.q_lwr_95 << ","
      << u.q_upr_95 << "\n";
}

} // namespace opakapaka_example

// Compatibility aliases for the current Opakapaka driver.
using opakapaka_example::compute_log_q_uncertainty_report;
using opakapaka_example::LogQUncertaintyReport;
using opakapaka_example::write_confidence_intervals_csv;
using opakapaka_example::write_correlation_matrix_csv;
using opakapaka_example::write_covariance_matrix_csv;
using opakapaka_example::write_standard_errors_csv;
using opakapaka_example::write_uncertainty_summary_csv;
