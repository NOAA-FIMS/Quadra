#pragma once

#include "../diagnostics/bigeye_functional_analysis_diagnostics.hpp"
#include "../reference_points/bigeye_reference_points.hpp"
#include "bigeye_fit_reports.hpp"

#include "../../../../../core/optimizer.hpp"

#include <string>
#include <vector>

namespace pifsc_bigeye_tuna {

struct BigeyeReportPaths {
  std::string summary =
      "examples/NMFS/pifsc_bigeye_tuna/level5_selectivity_contrast/outputs/"
      "bigeye_level5_fit_summary.csv";
  std::string trajectory =
      "examples/NMFS/pifsc_bigeye_tuna/level5_selectivity_contrast/outputs/"
      "bigeye_level5_fitted_trajectory.csv";
  std::string residual_diagnostics =
      "examples/NMFS/pifsc_bigeye_tuna/level5_selectivity_contrast/outputs/"
      "bigeye_level5_residual_diagnostics.csv";
  std::string selectivity =
      "examples/NMFS/pifsc_bigeye_tuna/level5_selectivity_contrast/outputs/"
      "bigeye_level5_selectivity_at_age.csv";
  std::string recruitment_deviations =
      "examples/NMFS/pifsc_bigeye_tuna/level5_selectivity_contrast/outputs/"
      "bigeye_level5_recruitment_deviations.csv";
  std::string objective_components =
      "examples/NMFS/pifsc_bigeye_tuna/level5_selectivity_contrast/outputs/"
      "bigeye_level5_objective_components.csv";
  std::string laplace_structure_text =
      "examples/NMFS/pifsc_bigeye_tuna/level5_selectivity_contrast/outputs/"
      "bigeye_level5_laplace_structure_report.txt";
  std::string laplace_structure_csv =
      "examples/NMFS/pifsc_bigeye_tuna/level5_selectivity_contrast/outputs/"
      "bigeye_level5_laplace_structure_report.csv";
  std::string reference_points =
      "examples/NMFS/pifsc_bigeye_tuna/level5_selectivity_contrast/outputs/"
      "bigeye_level5_reference_points.csv";
};

inline BigeyeReportPaths default_bigeye_report_paths() {
  return BigeyeReportPaths{};
}

template <class Objective>
inline void write_bigeye_report_suite(
    const BigeyeReportPaths &paths,
    const std::vector<Observation> &observations, Objective &objective,
    const quadra::ParameterVector &params, const quadra::OptResult &fit) {
  write_fit_summary(paths.summary, fit);
  write_fitted_trajectory(paths.trajectory, observations, fit);
  write_residual_diagnostics(paths.residual_diagnostics, observations, fit);
  write_selectivity_at_age(paths.selectivity, fit);
  write_recruitment_deviations(paths.recruitment_deviations, fit);
  write_objective_components(paths.objective_components, observations, fit);
  write_bigeye_reference_points_csv(paths.reference_points, observations, fit);
  write_bigeye_laplace_structure_report(paths.laplace_structure_text,
                                        paths.laplace_structure_csv, objective,
                                        params, fit);
}

} // namespace pifsc_bigeye_tuna

using pifsc_bigeye_tuna::BigeyeReportPaths;
using pifsc_bigeye_tuna::default_bigeye_report_paths;
using pifsc_bigeye_tuna::write_bigeye_report_suite;
