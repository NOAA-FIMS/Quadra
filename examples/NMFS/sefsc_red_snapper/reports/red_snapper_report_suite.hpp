#pragma once

#include "red_snapper_fit_reports.hpp"
#include "../diagnostics/red_snapper_functional_analysis_diagnostics.hpp"

#include "../../../../core/optimizer.hpp"

#include <string>
#include <vector>

namespace sefsc_red_snapper {

struct RedSnapperReportPaths
{
  std::string summary =
      "examples/NMFS/sefsc_red_snapper/outputs/quadra_fit_summary.csv";
  std::string trajectory =
      "examples/NMFS/sefsc_red_snapper/outputs/quadra_fitted_trajectory.csv";
  std::string residual_diagnostics =
      "examples/NMFS/sefsc_red_snapper/outputs/"
      "quadra_fit_residual_diagnostics.csv";
  std::string selectivity =
      "examples/NMFS/sefsc_red_snapper/outputs/selectivity_at_age.csv";
  std::string recruitment_deviations =
      "examples/NMFS/sefsc_red_snapper/outputs/recruitment_deviations.csv";
  std::string objective_components =
      "examples/NMFS/sefsc_red_snapper/outputs/"
      "quadra_fit_objective_components.csv";
  std::string laplace_structure_text =
      "examples/NMFS/sefsc_red_snapper/outputs/"
      "red_snapper_laplace_structure_report.txt";
  std::string laplace_structure_csv =
      "examples/NMFS/sefsc_red_snapper/outputs/"
      "red_snapper_laplace_structure_report.csv";
};

inline RedSnapperReportPaths default_red_snapper_report_paths()
{
  return RedSnapperReportPaths{};
}

template <class Objective>
inline void write_red_snapper_report_suite(
    const RedSnapperReportPaths &paths,
    const std::vector<Observation> &observations,
    Objective &objective,
    const quadra::ParameterVector &params,
    const quadra::OptResult &fit)
{
  write_fit_summary(paths.summary, fit);
  write_fitted_trajectory(paths.trajectory, observations, fit);
  write_residual_diagnostics(paths.residual_diagnostics, observations, fit);
  write_selectivity_at_age(paths.selectivity, fit);
  write_recruitment_deviations(paths.recruitment_deviations, fit);
  write_objective_components(paths.objective_components, observations, fit);
  write_red_snapper_laplace_structure_report(
      paths.laplace_structure_text,
      paths.laplace_structure_csv,
      objective, params, fit);
}

}  // namespace sefsc_red_snapper

using sefsc_red_snapper::RedSnapperReportPaths;
using sefsc_red_snapper::default_red_snapper_report_paths;
using sefsc_red_snapper::write_red_snapper_report_suite;
