#pragma once

#include "bigeye_safe_fixed_effect_diagnostics.hpp"
#include "../../../../../core/optimizer.hpp"

#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace pifsc_bigeye_tuna {

struct LonglineSlopeScanRow {
  int step_index = 0;
  double multiplier = std::numeric_limits<double>::quiet_NaN();
  double log_sel_slope_longline = std::numeric_limits<double>::quiet_NaN();
  double sel_slope_longline = std::numeric_limits<double>::quiet_NaN();
  double objective = std::numeric_limits<double>::quiet_NaN();
  double delta_objective = std::numeric_limits<double>::quiet_NaN();
  std::string status = "not_run";
  std::string message;
};

struct LonglineSlopeGeometryScan {
  std::size_t parameter_index = 5;
  std::string parameter_name = "log_sel_slope_longline";
  double base_log_slope = std::numeric_limits<double>::quiet_NaN();
  double base_slope = std::numeric_limits<double>::quiet_NaN();
  double base_objective = std::numeric_limits<double>::quiet_NaN();
  std::vector<LonglineSlopeScanRow> rows;
};

template <class Objective>
LonglineSlopeGeometryScan make_longline_slope_geometry_scan(
    Objective &objective, const quadra::ParameterVector &params,
    const quadra::OptResult &fit, quadra::LaplaceOptions opts,
    std::size_t parameter_index = 5) {
  LonglineSlopeGeometryScan scan;
  scan.parameter_index = parameter_index;

  if (fit.par.size() <= parameter_index) {
    LonglineSlopeScanRow row;
    row.status = "failed";
    row.message = "fit.par does not contain requested parameter index";
    scan.rows.push_back(row);
    return scan;
  }

  const auto names = safe_fixed_effect_names(params, fit);
  if (parameter_index < names.size()) scan.parameter_name = names[parameter_index];

  scan.base_log_slope = fit.par[parameter_index];
  scan.base_slope = std::exp(scan.base_log_slope);

  const auto base =
      safe_profiled_laplace_value_at_fixed(objective, params, fit.par, opts);
  scan.base_objective = base.value;

  const std::vector<double> multipliers = {
      0.40, 0.50, 0.60, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95,
      0.975, 1.00, 1.025, 1.05, 1.10, 1.15, 1.20, 1.30, 1.40,
      1.60, 1.80, 2.00};

  for (std::size_t i = 0; i < multipliers.size(); ++i) {
    LonglineSlopeScanRow row;
    row.step_index = static_cast<int>(i);
    row.multiplier = multipliers[i];
    row.sel_slope_longline = scan.base_slope * row.multiplier;
    row.log_sel_slope_longline = std::log(row.sel_slope_longline);

    auto theta = fit.par;
    theta[parameter_index] = row.log_sel_slope_longline;

    const auto val =
        safe_profiled_laplace_value_at_fixed(objective, params, theta, opts);

    row.objective = val.value;
    row.status = val.ok ? "ok" : "failed";
    row.message = val.message;

    if (val.ok && base.ok) row.delta_objective = val.value - base.value;

    scan.rows.push_back(row);
  }

  return scan;
}

inline void write_longline_slope_geometry_scan_csv(
    const LonglineSlopeGeometryScan &scan, const std::string &path) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("Could not open longline slope geometry scan CSV: " + path);

  out << "parameter_index,parameter_name,base_log_slope,base_slope,"
      << "base_objective,step_index,multiplier,log_sel_slope_longline,"
      << "sel_slope_longline,objective,delta_objective,status,message\n";

  out << std::setprecision(15);
  for (const auto &row : scan.rows) {
    out << scan.parameter_index << "," << scan.parameter_name << ","
        << scan.base_log_slope << "," << scan.base_slope << ","
        << scan.base_objective << "," << row.step_index << ","
        << row.multiplier << "," << row.log_sel_slope_longline << ","
        << row.sel_slope_longline << "," << row.objective << ","
        << row.delta_objective << "," << row.status << ","
        << '"' << row.message << '"' << "\n";
  }
}

inline void write_longline_slope_geometry_scan_text(
    const LonglineSlopeGeometryScan &scan, const std::string &path) {
  std::ofstream out(path);
  if (!out) throw std::runtime_error("Could not open longline slope geometry scan text: " + path);

  out << "Longline Selectivity Slope Geometry Scan\n";
  out << "=======================================\n\n";
  out << std::setprecision(15);
  out << "parameter_index:       " << scan.parameter_index << "\n";
  out << "parameter_name:        " << scan.parameter_name << "\n";
  out << "base_log_slope:        " << scan.base_log_slope << "\n";
  out << "base_slope:            " << scan.base_slope << "\n";
  out << "base_objective:        " << scan.base_objective << "\n\n";

  out << "Interpretation\n";
  out << "--------------\n";
  out << "This scan holds all fixed effects at the Level 6 fit except longline\n";
  out << "log-selectivity-slope. For each slope multiplier, random effects are\n";
  out << "re-solved and the profiled Laplace objective is evaluated. Failed rows\n";
  out << "mark local Huu/Laplace stability boundaries.\n\n";

  out << "Rows\n";
  out << "----\n";
  out << "step_index,multiplier,log_sel_slope_longline,sel_slope_longline,"
      << "objective,delta_objective,status,message\n";

  for (const auto &row : scan.rows) {
    out << row.step_index << "," << row.multiplier << ","
        << row.log_sel_slope_longline << "," << row.sel_slope_longline << ","
        << row.objective << "," << row.delta_objective << ","
        << row.status << "," << row.message << "\n";
  }
}

template <class Objective>
void write_longline_slope_geometry_scan(
    const std::string &text_path, const std::string &csv_path,
    Objective &objective, const quadra::ParameterVector &params,
    const quadra::OptResult &fit, quadra::LaplaceOptions opts,
    std::size_t parameter_index = 5) {
  const auto scan = make_longline_slope_geometry_scan(
      objective, params, fit, opts, parameter_index);

  write_longline_slope_geometry_scan_text(scan, text_path);
  write_longline_slope_geometry_scan_csv(scan, csv_path);
}

} // namespace pifsc_bigeye_tuna

using pifsc_bigeye_tuna::LonglineSlopeGeometryScan;
using pifsc_bigeye_tuna::write_longline_slope_geometry_scan;
