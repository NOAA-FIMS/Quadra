#pragma once

#include "../objective/bigeye_quadra_objective.hpp"
#include "../quadra/bigeye_age_structured.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <vector>

namespace pifsc_bigeye_tuna {

struct BigeyeFitReportPaths {
  std::string summary =
      "examples/NMFS/pifsc_bigeye_tuna/level18_longline_selectivity_diagnostic/"
      "outputs/bigeye_level18_fit_summary.csv";
  std::string trajectory =
      "examples/NMFS/pifsc_bigeye_tuna/level18_longline_selectivity_diagnostic/"
      "outputs/bigeye_level18_fitted_trajectory.csv";
  std::string residual_diagnostics =
      "examples/NMFS/pifsc_bigeye_tuna/level18_longline_selectivity_diagnostic/"
      "outputs/bigeye_level18_residual_diagnostics.csv";
  std::string selectivity =
      "examples/NMFS/pifsc_bigeye_tuna/level18_longline_selectivity_diagnostic/"
      "outputs/bigeye_level18_selectivity_at_age.csv";
  std::string recruitment_deviations =
      "examples/NMFS/pifsc_bigeye_tuna/level18_longline_selectivity_diagnostic/"
      "outputs/bigeye_level18_recruitment_deviations.csv";
  std::string objective_components =
      "examples/NMFS/pifsc_bigeye_tuna/level18_longline_selectivity_diagnostic/"
      "outputs/bigeye_level18_objective_components.csv";
  std::string laplace_structure_text =
      "examples/NMFS/pifsc_bigeye_tuna/level18_longline_selectivity_diagnostic/"
      "outputs/bigeye_level18_laplace_structure_report.txt";
  std::string laplace_structure_csv =
      "examples/NMFS/pifsc_bigeye_tuna/level18_longline_selectivity_diagnostic/"
      "outputs/bigeye_level18_laplace_structure_report.csv";
  std::string reference_points =
      "examples/NMFS/pifsc_bigeye_tuna/level18_longline_selectivity_diagnostic/"
      "outputs/bigeye_level18_reference_points.csv";
};

inline BigeyeFitReportPaths default_bigeye_report_paths() { return {}; }

inline double level11_direct_joint(const BigeyeQuadraObjective &objective,
                                   const quadra::OptResult &fit) {
  std::vector<double> full;
  for (double v : fit.par)
    full.push_back(v);
  for (double v : fit.u_hat)
    full.push_back(v);
  return const_cast<BigeyeQuadraObjective &>(objective)(full);
}

inline void write_bigeye_fit_summary(const std::string &path,
                                     const quadra::OptResult &fit) {
  std::ofstream out(path);
  if (!out)
    throw std::runtime_error("Could not open Level 18 summary: " + path);
  out << std::setprecision(12);
  out << "field,value\n";
  out << "objective," << fit.value << "\n";
  out << "grad_norm," << fit.grad_norm << "\n";
  out << "iterations," << fit.iterations << "\n";
  out << "converged," << (fit.converged ? "yes" : "no") << "\n";
  out << "message," << fit.message << "\n";
  out << "random_effects," << fit.u_hat.size() << "\n";
  if (fit.par.size() >= 5) {
    out << "log_r0," << fit.par[0] << "\n";
    out << "r0," << std::exp(fit.par[0]) << "\n";
    out << "log_m_fixed," << std::log(0.45) << "\n";
    out << "m_fixed," << 0.45 << "\n";
    out << "log_fbar," << fit.par[1] << "\n";
    out << "fbar," << std::exp(fit.par[1]) << "\n";
    out << "log_q_longline_fixed," << std::log(0.00005) << "\n";
    out << "q_longline_fixed," << 0.00005 << "\n";
    out << "log_q_purse_seine," << fit.par[2] << "\n";
    out << "q_purse_seine," << std::exp(fit.par[2]) << "\n";
    out << "logit_sel_a50_longline," << fit.par[3] << "\n";
    out << "sel_a50_longline," << 1.0 + 9.0 / (1.0 + std::exp(-fit.par[3]))
        << "\n";
    out << "log_sel_slope_longline," << fit.par[4] << "\n";
    out << "sel_slope_longline," << std::exp(fit.par[4]) << "\n";
    for (int a = 0; a < pifsc_bigeye_tuna::kAges; ++a) {
      const int idx = 5 + a;
      if (idx < static_cast<int>(fit.par.size())) {
        out << "init_log_number_dev_age_" << (a + 1) << "," << fit.par[idx]
            << "\n";
        out << "init_number_multiplier_age_" << (a + 1) << ","
            << std::exp(fit.par[idx]) << "\n";
      }
    }

    if (fit.par.size() >=
        static_cast<std::size_t>(5 + 2 * pifsc_bigeye_tuna::kAges)) {
      for (int a = 0; a < pifsc_bigeye_tuna::kAges; ++a) {
        const int idx = 5 + pifsc_bigeye_tuna::kAges + a;
        const double sel = 1.0 / (1.0 + std::exp(-fit.par[idx]));
        out << "logit_sel_purse_seine_age_" << (a + 1) << "," << fit.par[idx]
            << "\n";
        out << "sel_purse_seine_age_" << (a + 1) << "," << sel << "\n";
      }
    }
  }
}

inline void
write_bigeye_recruitment_deviations(const std::string &path,
                                    const BigeyeQuadraObjective &objective,
                                    const quadra::OptResult &fit) {
  std::ofstream out(path);
  if (!out)
    throw std::runtime_error("Could not open Level 18 recruitment: " + path);
  const auto years = objective.unique_years();
  out << std::setprecision(12);
  out << "year,rec_dev,recruitment_multiplier\n";
  for (std::size_t t = 0; t < years.size() && t < fit.u_hat.size(); ++t)
    out << years[t] << "," << fit.u_hat[t] << "," << std::exp(fit.u_hat[t])
        << "\n";
}

inline void
write_bigeye_objective_components(const std::string &path,
                                  const BigeyeQuadraObjective &objective,
                                  const quadra::OptResult &fit) {
  std::ofstream out(path);
  if (!out)
    throw std::runtime_error("Could not open Level 18 components: " + path);
  out << std::setprecision(12);
  out << "component,value\n";
  out << "joint_total," << level11_direct_joint(objective, fit) << "\n";
}

inline void write_bigeye_selectivity_at_age(const std::string &path,
                                            const quadra::OptResult &fit) {
  std::ofstream out(path);
  if (!out)
    throw std::runtime_error("Could not open Level 18 selectivity: " + path);
  out << "age,longline_selectivity,purse_seine_selectivity\n";
  if (fit.par.size() < 5)
    return;
  const double a50 = 1.0 + 9.0 / (1.0 + std::exp(-fit.par[3]));
  const double slope = std::exp(fit.par[4]);
  const double ps[10] = {1.00, 1.00, 0.85, 0.55, 0.25,
                         0.10, 0.04, 0.02, 0.01, 0.005};
  for (int a = 0; a < kAges; ++a) {
    const double ll =
        logistic_selectivity(static_cast<double>(a + 1), a50, slope);
    out << (a + 1) << "," << ll << "," << ps[a] << "\n";
  }
}

inline void
write_bigeye_fitted_trajectory(const std::string &path,
                               const BigeyeQuadraObjective &objective,
                               const quadra::OptResult &fit) {
  std::ofstream out(path);
  if (!out)
    throw std::runtime_error("Could not open Level 18 trajectory: " + path);
  const auto years = objective.unique_years();
  out << "year,rec_dev\n";
  for (std::size_t t = 0; t < years.size() && t < fit.u_hat.size(); ++t)
    out << years[t] << "," << fit.u_hat[t] << "\n";
}

inline void
write_bigeye_residual_diagnostics(const std::string &path,
                                  const BigeyeQuadraObjective &objective,
                                  const quadra::OptResult &fit) {
  std::ofstream out(path);
  if (!out)
    throw std::runtime_error("Could not open Level 18 residual diagnostics: " +
                             path);
  out << "metric,value,note\n";
  out << "direct_joint," << level11_direct_joint(objective, fit)
      << ",objective at fit.par+u_hat\n";
  out << "profiled_objective," << fit.value << ",optimizer fit value\n";
}

inline void write_bigeye_report_suite(const BigeyeFitReportPaths &paths,
                                      const std::vector<Observation> &,
                                      const BigeyeQuadraObjective &objective,
                                      const quadra::ParameterVector &,
                                      const quadra::OptResult &fit) {
  write_bigeye_fit_summary(paths.summary, fit);
  write_bigeye_fitted_trajectory(paths.trajectory, objective, fit);
  write_bigeye_residual_diagnostics(paths.residual_diagnostics, objective, fit);
  write_bigeye_selectivity_at_age(paths.selectivity, fit);
  write_bigeye_recruitment_deviations(paths.recruitment_deviations, objective,
                                      fit);
  write_bigeye_objective_components(paths.objective_components, objective, fit);
}

} // namespace pifsc_bigeye_tuna

using pifsc_bigeye_tuna::BigeyeFitReportPaths;
using pifsc_bigeye_tuna::default_bigeye_report_paths;
using pifsc_bigeye_tuna::write_bigeye_report_suite;
