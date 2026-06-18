#pragma once

#include "../quadra/red_snapper_age_structured.hpp"

#include "../../../../core/optimizer.hpp"

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>

namespace sefsc_red_snapper {

inline void write_fit_summary(const std::string &path, const quadra::OptResult &fit) {
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("Could not open fit summary CSV: " + path);
  }

  out << "field,value\n";
  out << std::setprecision(12);
  out << "objective," << fit.value << "\n";
  out << "joint_objective," << fit.joint_objective << "\n";
  out << "laplace_logdet," << fit.laplace_logdet << "\n";
  out << "laplace_constant," << fit.laplace_constant << "\n";
  out << "grad_norm," << fit.grad_norm << "\n";
  out << "iterations," << fit.iterations << "\n";
  out << "converged," << (fit.converged ? "yes" : "no") << "\n";
  out << "message," << fit.message << "\n";
  out << "laplace,yes\n";
  out << "random_effects," << fit.u_hat.size() << "\n";

  if (fit.par.size() >= 3) {
    out << "log_r0," << fit.par[0] << "\n";
    out << "r0," << std::exp(fit.par[0]) << "\n";
    out << "log_fbar," << fit.par[1] << "\n";
    out << "fbar," << std::exp(fit.par[1]) << "\n";
    out << "log_q," << fit.par[2] << "\n";
    out << "q," << std::exp(fit.par[2]) << "\n";
    if (fit.par.size() >= 5) {
      const double sel_a50 = 1.0 + 9.0 / (1.0 + std::exp(-fit.par[3]));
      const double sel_slope = std::exp(fit.par[4]);
      out << "logit_sel_a50," << fit.par[3] << "\n";
      out << "sel_a50," << sel_a50 << "\n";
      out << "log_sel_slope," << fit.par[4] << "\n";
      out << "sel_slope," << sel_slope << "\n";
    }
  }
}


inline void write_fitted_trajectory(
    const std::string &path,
    const std::vector<sefsc_red_snapper::Observation> &observations,
    const quadra::OptResult &fit) {
  if (fit.par.size() < 3) {
    throw std::runtime_error(
        "Cannot write fitted trajectory: expected at least 3 fixed parameters");
  }

  sefsc_red_snapper::AgeStructuredParams params;
  params.log_r0 = fit.par[0];
  params.log_fbar = fit.par[1];
  params.log_q = fit.par[2];
  if (fit.par.size() >= 5) {
    params.sel_a50 = 1.0 + 9.0 / (1.0 + std::exp(-fit.par[3]));
    params.sel_slope = std::exp(fit.par[4]);
  }

  const auto rows = sefsc_red_snapper::run_deterministic_age_structured_model(
      observations, params);

  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("Could not open fitted trajectory CSV: " + path);
  }

  out << "year,recruitment,total_biomass,ssb_proxy,depletion,Fbar,"
      << "catch_obs,catch_hat,catch_log_residual,index_obs,index_hat,"
      << "index_log_residual\n";

  out << std::fixed << std::setprecision(6);

  for (const auto &row : rows) {
    const double catch_log_residual =
        std::log(std::max(row.catch_obs, 1.0e-12)) -
        std::log(std::max(row.catch_hat, 1.0e-12));
    const double index_log_residual =
        std::log(std::max(row.index_obs, 1.0e-12)) -
        std::log(std::max(row.index_hat, 1.0e-12));

    out << row.year << "," << row.recruitment << "," << row.total_biomass << ","
        << row.ssb_proxy << "," << row.depletion << "," << row.fbar << ","
        << row.catch_obs << "," << row.catch_hat << "," << catch_log_residual
        << "," << row.index_obs << "," << row.index_hat << ","
        << index_log_residual << "\n";
  }
}

struct ResidualDiagnostics {
  int n = 0;
  double catch_rmse_log = 0.0;
  double index_rmse_log = 0.0;
  double catch_mean_log_residual = 0.0;
  double index_mean_log_residual = 0.0;
  double max_abs_catch_log_residual = 0.0;
  double max_abs_index_log_residual = 0.0;
};

inline void write_residual_diagnostics(
    const std::string &path,
    const std::vector<sefsc_red_snapper::Observation> &observations,
    const quadra::OptResult &fit) {
  sefsc_red_snapper::AgeStructuredParams params;
  params.log_r0 = fit.par[0];
  params.log_fbar = fit.par[1];
  params.log_q = fit.par[2];
  if (fit.par.size() >= 5) {
    params.sel_a50 = 1.0 + 9.0 / (1.0 + std::exp(-fit.par[3]));
    params.sel_slope = std::exp(fit.par[4]);
  }

  const auto rows = sefsc_red_snapper::run_deterministic_age_structured_model(
      observations, params);

  ResidualDiagnostics d;
  d.n = static_cast<int>(rows.size());

  double catch_sum = 0.0, catch_ss = 0.0;
  double index_sum = 0.0, index_ss = 0.0;

  for (const auto &row : rows) {
    const double cr = std::log(std::max(row.catch_obs, 1.0e-12)) -
                      std::log(std::max(row.catch_hat, 1.0e-12));
    const double ir = std::log(std::max(row.index_obs, 1.0e-12)) -
                      std::log(std::max(row.index_hat, 1.0e-12));

    catch_sum += cr;
    catch_ss += cr * cr;
    index_sum += ir;
    index_ss += ir * ir;

    d.max_abs_catch_log_residual =
        std::max(d.max_abs_catch_log_residual, std::abs(cr));
    d.max_abs_index_log_residual =
        std::max(d.max_abs_index_log_residual, std::abs(ir));
  }

  if (d.n > 0) {
    d.catch_mean_log_residual = catch_sum / d.n;
    d.index_mean_log_residual = index_sum / d.n;
    d.catch_rmse_log = std::sqrt(catch_ss / d.n);
    d.index_rmse_log = std::sqrt(index_ss / d.n);
  }

  std::ofstream out(path);
  out << "metric,value,note\n";
  out << std::setprecision(12);
  out << "n," << d.n << ",number of fitted years\n";
  out << "catch_rmse_log," << d.catch_rmse_log
      << ",root mean squared log catch residual\n";
  out << "index_rmse_log," << d.index_rmse_log
      << ",root mean squared log index residual\n";
  out << "catch_mean_log_residual," << d.catch_mean_log_residual
      << ",mean log observed minus predicted catch\n";
  out << "index_mean_log_residual," << d.index_mean_log_residual
      << ",mean log observed minus predicted index\n";
  out << "max_abs_catch_log_residual," << d.max_abs_catch_log_residual
      << ",maximum absolute log catch residual\n";
  out << "max_abs_index_log_residual," << d.max_abs_index_log_residual
      << ",maximum absolute log index residual\n";
}

inline void write_selectivity_at_age(const std::string &path,
                              const quadra::OptResult &fit) {
  if (fit.par.size() < 5) {
    return;
  }

  const double a50 = 1.0 + 9.0 / (1.0 + std::exp(-fit.par[3]));
  const double slope = std::exp(fit.par[4]);

  std::ofstream out(path);
  out << "age,selectivity\n";

  for (int age = 1; age <= sefsc_red_snapper::kAges; ++age) {
    const double sel = 1.0 / (1.0 + std::exp(-slope * (age - a50)));
    out << age << "," << sel << "\n";
  }
}

inline void write_recruitment_deviations(const std::string &path,
                                  const quadra::OptResult &fit) {
  std::ofstream out(path);
  out << "year,log_rec_dev,rec_multiplier\n";
  out << std::setprecision(12);

  for (std::size_t i = 0; i < fit.u_hat.size(); ++i) {
    const double u = fit.u_hat[i];
    out << (i + 1) << "," << u << "," << std::exp(u) << "\n";
  }
}

inline void write_objective_components(
    const std::string &path,
    const std::vector<sefsc_red_snapper::Observation> &observations,
    const quadra::OptResult &fit) {
  if (fit.par.size() < 5 || fit.u_hat.size() < observations.size()) {
    throw std::runtime_error(
        "Cannot write objective components: missing fit values");
  }

  const double log_r0 = fit.par[0];
  const double log_fbar = fit.par[1];
  const double log_q = fit.par[2];
  const double logit_sel_a50 = fit.par[3];
  const double log_sel_slope = fit.par[4];

  const double r0 = std::exp(log_r0);
  const double m = 0.18;
  const double fbar = std::exp(log_fbar);
  const double q = std::exp(log_q);
  const double sel_a50 = 1.0 + 9.0 / (1.0 + std::exp(-logit_sel_a50));
  const double sel_slope = std::exp(log_sel_slope);

  const double sigma_log_index = 0.20;
  const double sigma_log_catch = 0.15;
  const double sigma_rec_dev = 0.35;
  const double age_comp_effective_n = 2.0;
  const double min_positive = 1.0e-12;

  const auto weight = sefsc_red_snapper::default_weight_at_age();

  std::array<double, sefsc_red_snapper::kAges> selectivity{};
  for (int a = 0; a < sefsc_red_snapper::kAges; ++a) {
    selectivity[static_cast<std::size_t>(a)] =
        sefsc_red_snapper::logistic_selectivity(static_cast<double>(a + 1),
                                                sel_a50, sel_slope);
  }

  std::array<double, sefsc_red_snapper::kAges> n{};
  n[0] = r0;
  for (int a = 1; a < sefsc_red_snapper::kAges; ++a) {
    n[static_cast<std::size_t>(a)] =
        n[static_cast<std::size_t>(a - 1)] * std::exp(-m);
  }
  n[static_cast<std::size_t>(sefsc_red_snapper::kAges - 1)] =
      n[static_cast<std::size_t>(sefsc_red_snapper::kAges - 1)] /
      (1.0 - std::exp(-m));

  auto normal_prior = [](double x, double mean, double sd) {
    const double z = (x - mean) / sd;
    return 0.5 * z * z;
  };

  double fixed_prior_nll = 0.0;
  double rec_prior_nll = 0.0;
  double index_nll = 0.0;
  double catch_nll = 0.0;
  double age_comp_nll = 0.0;

  fixed_prior_nll += normal_prior(log_r0, std::log(1200.0), 1.0);
  fixed_prior_nll += normal_prior(log_fbar, std::log(0.025), 0.75);
  fixed_prior_nll += normal_prior(log_q, std::log(0.00005), 1.0);
  fixed_prior_nll += normal_prior(sel_a50, 4.0, 0.75);
  fixed_prior_nll += normal_prior(log_sel_slope, std::log(1.2), 0.35);

  for (std::size_t t = 0; t < observations.size(); ++t) {
    const auto &obs = observations[t];
    const double rec_dev = fit.u_hat[t];

    rec_prior_nll += 0.5 * std::pow(rec_dev / sigma_rec_dev, 2.0);

    double biomass = 0.0;
    for (int a = 0; a < sefsc_red_snapper::kAges; ++a) {
      biomass +=
          n[static_cast<std::size_t>(a)] * weight[static_cast<std::size_t>(a)];
    }

    double catch_hat = 0.0;
    for (int a = 0; a < sefsc_red_snapper::kAges; ++a) {
      const auto i = static_cast<std::size_t>(a);
      const double f_a = fbar * selectivity[i];
      const double z_a = m + f_a;
      const double harvest_rate = (f_a / z_a) * (1.0 - std::exp(-z_a));
      catch_hat += n[i] * weight[i] * harvest_rate;
    }

    const double index_hat = q * biomass;

    if (obs.index > 0.0) {
      const double z =
          (std::log(obs.index) - std::log(std::max(index_hat, min_positive))) /
          sigma_log_index;
      index_nll += 0.5 * z * z;
    }

    if (obs.catch_mt > 0.0) {
      const double z = (std::log(obs.catch_mt) -
                        std::log(std::max(catch_hat, min_positive))) /
                       sigma_log_catch;
      catch_nll += 0.5 * z * z;
    }

    std::array<double, sefsc_red_snapper::kAges> pred_age_comp{};
    double selected_numbers_sum = 0.0;
    for (int a = 0; a < sefsc_red_snapper::kAges; ++a) {
      const auto i = static_cast<std::size_t>(a);
      pred_age_comp[i] = n[i] * selectivity[i];
      selected_numbers_sum += pred_age_comp[i];
    }

    for (int a = 0; a < sefsc_red_snapper::kAges; ++a) {
      const auto i = static_cast<std::size_t>(a);
      pred_age_comp[i] =
          pred_age_comp[i] / std::max(selected_numbers_sum, min_positive);

      const double obs_a = std::max(obs.age_comp[i], 0.0);
      if (obs_a > 0.0) {
        age_comp_nll -= age_comp_effective_n * obs_a *
                        std::log(std::max(pred_age_comp[i], min_positive));
      }
    }

    std::array<double, sefsc_red_snapper::kAges> next{};
    next[0] = r0 * std::exp(rec_dev);
    for (int a = 1; a < sefsc_red_snapper::kAges; ++a) {
      const auto prev = static_cast<std::size_t>(a - 1);
      const auto cur = static_cast<std::size_t>(a);
      const double f_prev = fbar * selectivity[prev];
      const double z_prev = m + f_prev;
      next[cur] = n[prev] * std::exp(-z_prev);
    }

    const int plus_group = sefsc_red_snapper::kAges - 1;
    const auto pg = static_cast<std::size_t>(plus_group);
    const double f_pg = fbar * selectivity[pg];
    const double z_pg = m + f_pg;
    next[pg] += n[pg] * std::exp(-z_pg);

    n = next;
  }

  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("Could not open component CSV: " + path);
  }

  out << "component,value\n";
  out << std::setprecision(12);
  out << "fixed_prior_nll," << fixed_prior_nll << "\n";
  out << "rec_prior_nll," << rec_prior_nll << "\n";
  out << "index_nll," << index_nll << "\n";
  out << "catch_nll," << catch_nll << "\n";
  out << "age_comp_nll," << age_comp_nll << "\n";
  out << "joint_total,"
      << fixed_prior_nll + rec_prior_nll + index_nll + catch_nll + age_comp_nll
      << "\n";
}

}  // namespace sefsc_red_snapper

using sefsc_red_snapper::write_fit_summary;
using sefsc_red_snapper::write_fitted_trajectory;
using sefsc_red_snapper::write_objective_components;
using sefsc_red_snapper::write_recruitment_deviations;
using sefsc_red_snapper::write_residual_diagnostics;
using sefsc_red_snapper::write_selectivity_at_age;
