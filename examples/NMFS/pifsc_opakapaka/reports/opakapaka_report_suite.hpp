#pragma once

#include "../diagnostics/opakapaka_biomass_covariance_diagnostics.hpp"
#include "../diagnostics/opakapaka_logq_diagnostics.hpp"
#include "../diagnostics/opakapaka_projection_uncertainty.hpp"
#include "../diagnostics/opakapaka_random_effect_diagnostics.hpp"
#include "../quadra/opakapaka_model.hpp"

#include "../../../../core/uncertainty/reporting.hpp"
#include "../../../../core/uncertainty/selected_inverse_diagonal.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>
#include <utility>
#include <vector>

namespace opakapaka_example {

inline void write_derived_quantities_csv(
  const std::string &path,
  const std::vector<opakapaka_example::Observation> &data,
  const std::vector<double> &u_hat, double q_hat)
{
std::ofstream out(path);
out << "year,biomass,index_hat,depletion,F_proxy\n";
const double b0 = u_hat.empty() ? std::numeric_limits<double>::quiet_NaN()
                                : std::exp(u_hat.front());
for (std::size_t i = 0; i < data.size() && i < u_hat.size(); ++i)
{
  const double biomass = std::exp(u_hat[i]);
  const double depletion =
      b0 > 0.0 ? biomass / b0 : std::numeric_limits<double>::quiet_NaN();
  const double f_proxy = biomass > 0.0
                             ? data[i].catch_mt / biomass
                             : std::numeric_limits<double>::quiet_NaN();
  out << data[i].year << "," << biomass << "," << q_hat * biomass << ","
      << depletion << "," << f_proxy << "\n";
}
}


inline void write_derived_quantity_uncertainty_csv(
  const std::string &path,
  const std::vector<opakapaka_example::Observation> &data,
  const std::vector<double> &u_hat, double q_hat,
  const quadra::uncertainty::SelectedInverseDiagonalResult &u_cov,
  const Eigen::SparseMatrix<double> &h_uu)
{
std::ofstream out(path);
out << "year,quantity,estimate,se,lwr_95,upr_95,note\n";

if (u_hat.empty() || data.empty())
{
  return;
}

const double b0 = std::exp(u_hat.front());
const double var_log_b0 = (u_cov.success && !u_cov.variance.empty())
                              ? u_cov.variance.front()
                              : std::numeric_limits<double>::quiet_NaN();

// QUADRA_OPAKAPAKA_DEPLETION_COVARIANCE_PAIRS_V1
// Request Cov(log_B[t], log_B[0]) so depletion uncertainty uses:
// Var(log(B_t/B_0)) = Var(log_B_t) + Var(log_B_0) - 2 Cov(log_B_t, log_B_0).
std::vector<std::pair<int, int>> depletion_covariance_pairs;
depletion_covariance_pairs.reserve(u_hat.size());
for (std::size_t i = 0; i < u_hat.size(); ++i)
{
  depletion_covariance_pairs.emplace_back(static_cast<int>(i), 0);
}

const auto depletion_covariances =
    quadra::uncertainty::selected_inverse_entries_from_spd_hessian(
        h_uu, depletion_covariance_pairs);

for (std::size_t i = 0; i < data.size() && i < u_hat.size(); ++i)
{
  const double log_b = u_hat[i];
  const double biomass = std::exp(log_b);
  const double index_hat = q_hat * biomass;
  const double depletion =
      b0 > 0.0 ? biomass / b0 : std::numeric_limits<double>::quiet_NaN();
  const double f_proxy = biomass > 0.0
                             ? data[i].catch_mt / biomass
                             : std::numeric_limits<double>::quiet_NaN();

  const double var_log_b = (u_cov.success && i < u_cov.variance.size())
                               ? u_cov.variance[i]
                               : std::numeric_limits<double>::quiet_NaN();

  const double se_biomass = (std::isfinite(var_log_b) && var_log_b >= 0.0)
                                ? biomass * std::sqrt(var_log_b)
                                : std::numeric_limits<double>::quiet_NaN();

  const double se_index = (std::isfinite(var_log_b) && var_log_b >= 0.0)
                              ? index_hat * std::sqrt(var_log_b)
                              : std::numeric_limits<double>::quiet_NaN();

  double cov_log_b_i_b0 = std::numeric_limits<double>::quiet_NaN();
  if (depletion_covariances.success &&
      i < depletion_covariances.entries.size())
  {
    cov_log_b_i_b0 = depletion_covariances.entries[i].covariance;
  }

  const double var_log_depletion =
      (std::isfinite(var_log_b) && std::isfinite(var_log_b0) &&
       std::isfinite(cov_log_b_i_b0))
          ? var_log_b + var_log_b0 - 2.0 * cov_log_b_i_b0
          : std::numeric_limits<double>::quiet_NaN();

  const double se_depletion =
      (std::isfinite(var_log_depletion) && var_log_depletion >= 0.0)
          ? depletion * std::sqrt(var_log_depletion)
          : std::numeric_limits<double>::quiet_NaN();

  const double se_f_proxy = (std::isfinite(var_log_b) && var_log_b >= 0.0)
                                ? f_proxy * std::sqrt(var_log_b)
                                : std::numeric_limits<double>::quiet_NaN();

  auto write_row = [&](const char *quantity, double estimate, double se,
                       const char *note)
  {
    const double lwr = std::isfinite(se)
                           ? estimate - 1.96 * se
                           : std::numeric_limits<double>::quiet_NaN();
    const double upr = std::isfinite(se)
                           ? estimate + 1.96 * se
                           : std::numeric_limits<double>::quiet_NaN();
    out << data[i].year << "," << quantity << "," << estimate << "," << se
        << "," << lwr << "," << upr << "," << note << "\n";
  };

  write_row("biomass", biomass, se_biomass,
            "level1_delta_method_conditional_random_effect_diagonal");
  write_row("index_hat", index_hat, se_index,
            "level1_delta_method_conditional_random_effect_diagonal");
  write_row("depletion", depletion, se_depletion,
            "level1_delta_method_selected_inverse_cov_logBt_logB0");
  write_row("F_proxy", f_proxy, se_f_proxy,
            "level1_delta_method_conditional_random_effect_diagonal");
}
}


inline void write_derived_quantity_correlation_csv(
  const std::string &path,
  const std::vector<opakapaka_example::Observation> &data,
  const quadra::uncertainty::SelectedInverseDiagonalResult &u_cov,
  const quadra::uncertainty::SelectedInverseEntriesResult
      &depletion_covariances)
{
std::ofstream out(path);
out << "year,variance_logB0,variance_logBt,covariance_logBt_logB0,"
    << "correlation_logBt_logB0,note\n";

const double var_log_b0 = (u_cov.success && !u_cov.variance.empty())
                              ? u_cov.variance.front()
                              : std::numeric_limits<double>::quiet_NaN();

const std::size_t n = std::min(data.size(), u_cov.variance.size());

for (std::size_t i = 0; i < n; ++i)
{
  const double var_log_bt = u_cov.variance[i];

  double cov_log_bt_b0 = std::numeric_limits<double>::quiet_NaN();
  if (depletion_covariances.success &&
      i < depletion_covariances.entries.size())
  {
    cov_log_bt_b0 = depletion_covariances.entries[i].covariance;
  }

  double corr = std::numeric_limits<double>::quiet_NaN();
  if (std::isfinite(var_log_b0) && std::isfinite(var_log_bt) &&
      std::isfinite(cov_log_bt_b0) && var_log_b0 > 0.0 && var_log_bt > 0.0)
  {
    corr = cov_log_bt_b0 / std::sqrt(var_log_b0 * var_log_bt);

    // Guard tiny numerical drift outside [-1, 1].
    if (corr > 1.0 && corr < 1.0 + 1.0e-10)
      corr = 1.0;
    if (corr < -1.0 && corr > -1.0 - 1.0e-10)
      corr = -1.0;
  }

  out << data[i].year << "," << var_log_b0 << "," << var_log_bt << ","
      << cov_log_bt_b0 << "," << corr << ","
      << "selected_inverse_covariance_diagnostic_logBt_logB0\n";
}
}


inline void write_runtime_memory_summary_csv(const std::string &path,
                                           double runtime_ms,
                                           std::size_t random_effects,
                                           std::size_t hessian_nonzeros)
{
std::ofstream out(path);
out << "field,value\n";
out << "fit_runtime_ms," << runtime_ms << "\n";
out << "random_effects," << random_effects << "\n";
out << "hessian_nonzeros," << hessian_nonzeros << "\n";
out << "peak_rss_mb,\n";
out << "note,peak RSS is captured by benchmark runner rather than model "
       "executable\n";
}



template <class Model>
inline void write_opakapaka_report_suite(
    Model &model,
    quadra::ParameterVector &params,
    quadra::LaplaceOptions &opts,
    const quadra::OptResult &fit,
    const std::vector<Observation> &data,
    const std::vector<ProjectionRow> &projection,
    const Eigen::SparseMatrix<double> &final_h_uu)
{
  write_fit_summary_csv(
      "examples/NMFS/pifsc_opakapaka/outputs/synthetic_fit_summary.csv", fit);

  const auto logq_uncertainty =
      compute_log_q_uncertainty_report(model, params, opts, fit);

  write_uncertainty_summary_csv(
      "examples/NMFS/pifsc_opakapaka/outputs/uncertainty_summary.csv",
      logq_uncertainty);
  write_covariance_matrix_csv(
      "examples/NMFS/pifsc_opakapaka/outputs/covariance_matrix.csv",
      logq_uncertainty);
  write_correlation_matrix_csv(
      "examples/NMFS/pifsc_opakapaka/outputs/correlation_matrix.csv");
  write_standard_errors_csv(
      "examples/NMFS/pifsc_opakapaka/outputs/standard_errors.csv",
      logq_uncertainty);
  write_confidence_intervals_csv(
      "examples/NMFS/pifsc_opakapaka/outputs/confidence_intervals.csv",
      logq_uncertainty);

  write_random_effect_uncertainty_csv(
      "examples/NMFS/pifsc_opakapaka/outputs/random_effect_uncertainty.csv",
      fit.u_hat, final_h_uu);

  write_derived_quantities_csv(
      "examples/NMFS/pifsc_opakapaka/outputs/derived_quantities.csv", data,
      fit.u_hat, std::exp(fit.par.at(0)));

  const auto random_effect_covariance_diag =
      quadra::uncertainty::selected_inverse_diagonal_from_spd_hessian(
          final_h_uu);

  write_derived_quantity_uncertainty_csv(
      "examples/NMFS/pifsc_opakapaka/outputs/derived_quantity_uncertainty.csv",
      data, fit.u_hat, std::exp(fit.par.at(0)), random_effect_covariance_diag,
      final_h_uu);

  {
    std::vector<std::pair<int, int>> depletion_covariance_pairs;
    depletion_covariance_pairs.reserve(fit.u_hat.size());
    for (std::size_t i = 0; i < fit.u_hat.size(); ++i)
    {
      depletion_covariance_pairs.emplace_back(static_cast<int>(i), 0);
    }

    const auto depletion_covariances =
        quadra::uncertainty::selected_inverse_entries_from_spd_hessian(
            final_h_uu, depletion_covariance_pairs);

    write_derived_quantity_correlation_csv(
        "examples/NMFS/pifsc_opakapaka/outputs/"
        "derived_quantity_correlation.csv",
        data, random_effect_covariance_diag, depletion_covariances);
  }

  write_biomass_covariance_matrix_csv(
      "examples/NMFS/pifsc_opakapaka/outputs/biomass_covariance_matrix.csv",
      data, fit.u_hat, final_h_uu);

  write_biomass_correlation_matrix_csv(
      "examples/NMFS/pifsc_opakapaka/outputs/biomass_correlation_matrix.csv",
      data, fit.u_hat, final_h_uu);

  write_biomass_covariance_diagnostics_csv(
      "examples/NMFS/pifsc_opakapaka/outputs/"
      "biomass_covariance_diagnostics.csv",
      data, fit.u_hat, final_h_uu);

  write_biomass_correlation_decay_csv(
      "examples/NMFS/pifsc_opakapaka/outputs/biomass_correlation_decay.csv",
      data, fit.u_hat, final_h_uu);

  {
    const std::size_t n = std::min(data.size(), fit.u_hat.size());
    const Eigen::MatrixXd log_b_cov =
        compute_log_b_covariance_submatrix(data, fit.u_hat, final_h_uu);
    Eigen::VectorXd log_b_core(static_cast<Eigen::Index>(n));
    for (std::size_t i = 0; i < n; ++i)
    {
      log_b_core[static_cast<Eigen::Index>(i)] = fit.u_hat[i];
    }
    const Eigen::MatrixXd biomass_cov_core =
        quadra::uncertainty::lognormal_delta_covariance(log_b_core,
                                                        log_b_cov);
    const Eigen::MatrixXd biomass_corr_core =
        quadra::uncertainty::covariance_to_correlation_matrix(biomass_cov_core);
    const auto biomass_diag_core =
        quadra::uncertainty::diagnose_covariance_matrix(biomass_cov_core);
    quadra::uncertainty::write_covariance_diagnostics_csv(
        "examples/NMFS/pifsc_opakapaka/outputs/"
        "biomass_covariance_diagnostics_core.csv",
        biomass_diag_core);
    const auto biomass_decay_core =
        quadra::uncertainty::correlation_decay_summary(biomass_corr_core);
    quadra::uncertainty::write_correlation_decay_csv(
        "examples/NMFS/pifsc_opakapaka/outputs/"
        "biomass_correlation_decay_core.csv",
        biomass_decay_core);
  }

  const double terminal_log_b_variance =
      (!random_effect_covariance_diag.variance.empty())
          ? random_effect_covariance_diag.variance.back()
          : std::numeric_limits<double>::quiet_NaN();

  write_projection_uncertainty_envelopes_csv(
      "examples/NMFS/pifsc_opakapaka/outputs/projection_uncertainty.csv",
      projection, fit.u_hat, std::exp(fit.par.at(0)), terminal_log_b_variance,
      1000);

  write_runtime_memory_summary_csv(
      "examples/NMFS/pifsc_opakapaka/outputs/runtime_memory_summary.csv",
      std::numeric_limits<double>::quiet_NaN(), fit.u_hat.size(), 58);

  write_projection_csv("examples/NMFS/pifsc_opakapaka/outputs/"
                       "synthetic_projection_scenarios.csv",
                       projection);
}

}  // namespace opakapaka_example

using opakapaka_example::write_opakapaka_report_suite;
