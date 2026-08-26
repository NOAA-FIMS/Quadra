#ifndef QUADRA_TUNA_ASSESSMENT_ACCEPTANCE_HPP
#define QUADRA_TUNA_ASSESSMENT_ACCEPTANCE_HPP
#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "core/laplace/joint_only_exact_lbfgs_optimizer.hpp"
#include "core/laplace/laplace_exact_lbfgs_optimizer.hpp"
#include "core/model/model_context.hpp"
#include "core/model/parameter_partition.hpp"
#include "include/quadra/stats/laplace.hpp"
#include "tuna_spatial_assessment_model.hpp"

namespace quadra {

struct TunaAssessmentRunSummary {
  std::string label_m;
  int n_years_m = 0;
  double nll_m = std::numeric_limits<double>::quiet_NaN();
  double r0_m = std::numeric_limits<double>::quiet_NaN();
  double steepness_m = std::numeric_limits<double>::quiet_NaN();
  double ssb0_m = std::numeric_limits<double>::quiet_NaN();
  double ssb_terminal_m = std::numeric_limits<double>::quiet_NaN();
  double depletion_terminal_m = std::numeric_limits<double>::quiet_NaN();
  std::vector<double> ssb_by_year_m;
};

struct TunaObservationDiagnostic {
  std::string component_m;
  int year_m = 0;
  int season_m = 0;
  int fleet_m = 0;
  int region_m = 0;
  int age_m = 0;
  double observed_m = std::numeric_limits<double>::quiet_NaN();
  double predicted_m = std::numeric_limits<double>::quiet_NaN();
  double standard_deviation_m = std::numeric_limits<double>::quiet_NaN();
  double standardized_residual_m = std::numeric_limits<double>::quiet_NaN();
};

struct TunaResidualSummary {
  std::string component_m;
  int fleet_m = 0;
  int n_m = 0;
  double mean_residual_m = std::numeric_limits<double>::quiet_NaN();
  double sdnr_m = std::numeric_limits<double>::quiet_NaN();
  double rmse_m = std::numeric_limits<double>::quiet_NaN();
};

struct TunaStratifiedResidualSummary {
  std::string component_m;
  int fleet_m = 0;
  std::string stratum_m;
  int level_m = 0;
  int n_m = 0;
  double mean_residual_m = std::numeric_limits<double>::quiet_NaN();
  double sdnr_m = std::numeric_limits<double>::quiet_NaN();
  double rmse_m = std::numeric_limits<double>::quiet_NaN();
};

struct TunaLikelihoodComponent {
  std::string component_m;
  double nll_m = std::numeric_limits<double>::quiet_NaN();
  double share_absolute_m = std::numeric_limits<double>::quiet_NaN();
  int observation_units_m = 0;
  double nll_per_unit_m = std::numeric_limits<double>::quiet_NaN();
};

struct TunaParameterDiagnostic {
  std::string name_m;
  double value_m = std::numeric_limits<double>::quiet_NaN();
  double initial_m = std::numeric_limits<double>::quiet_NaN();
  double displacement_m = std::numeric_limits<double>::quiet_NaN();
  bool random_effect_m = false;
  bool extreme_m = false;
};

struct TunaCatchabilityAvailabilityDiagnostic {
  int fleet_m = 0;
  int season_m = 0;
  int region_m = 0;
  double log_fishing_q_m = std::numeric_limits<double>::quiet_NaN();
  double log_index_q_m = std::numeric_limits<double>::quiet_NaN();
  double log_availability_m = std::numeric_limits<double>::quiet_NaN();
  double log_fishing_q_availability_m =
      std::numeric_limits<double>::quiet_NaN();
};

struct TunaDiagnosticBundle {
  std::vector<TunaObservationDiagnostic> observations_m;
  std::vector<TunaResidualSummary> residual_summaries_m;
  std::vector<TunaStratifiedResidualSummary> stratified_residuals_m;
  std::vector<TunaLikelihoodComponent> likelihood_components_m;
  std::vector<TunaParameterDiagnostic> parameters_m;
  std::vector<TunaCatchabilityAvailabilityDiagnostic>
      catchability_availability_m;
};

struct TunaRetrospectivePoint {
  int peel_m = 0;
  TunaAssessmentRunSummary summary_m;
};

struct TunaRetrospectiveResult {
  std::vector<TunaRetrospectivePoint> points_m;
  double mohns_rho_ssb_m = std::numeric_limits<double>::quiet_NaN();
};

struct TunaSensitivityScenario {
  std::string label_m;
  double natural_mortality_multiplier_m = 1.0;
  double effort_multiplier_m = 1.0;
  double availability_multiplier_m = 1.0;
  double index_multiplier_m = 1.0;
};

struct TunaFitOptions {
  int multistart_m = 3;
  int max_iterations_per_phase_m = 50;
  int lbfgs_memory_m = 7;
  int hdot_workers_m = 0;
  int lbfgs_print_every_m = 0;
  int max_evaluations_per_phase_m = 0;
  double gradient_tolerance_m = 1e-4;
  double line_search_acceptance_gradient_tolerance_m = 5e-3;
  double initial_step_m = 1.0;
  double min_step_m = 1e-8;
  double jitter_sd_m = 0.35;
  std::uint64_t seed_m = 42;
  std::vector<std::string> fixed_parameter_names_m;
  std::vector<double> fixed_parameter_values_m;
  std::vector<TunaAssessmentPhase> phase_sequence_m = {
      TunaAssessmentPhase::InitializeRecruitment,
      TunaAssessmentPhase::InitializeCatchability,
      TunaAssessmentPhase::InitializeMovement, TunaAssessmentPhase::Full};
};

struct TunaFitResult {
  bool converged_m = false;
  int best_start_index_m = -1;
  int objective_evaluations_m = 0;
  int total_iterations_m = 0;
  double gradient_norm_m = std::numeric_limits<double>::quiet_NaN();
  std::string message_m;
  std::vector<double> best_parameters_m;
  TunaAssessmentRunSummary summary_m;
  struct PhaseDiagnostic {
    TunaAssessmentPhase phase_m = TunaAssessmentPhase::Full;
    int start_index_m = -1;
    int free_fixed_count_m = 0;
    int locked_fixed_count_m = 0;
    int iterations_m = 0;
    double laplace_objective_m = std::numeric_limits<double>::quiet_NaN();
    double gradient_norm_m = std::numeric_limits<double>::quiet_NaN();
    std::vector<std::string> largest_gradient_parameters_m;
    std::vector<double> largest_gradient_values_m;
    bool converged_m = false;
    bool logdet_ok_m = false;
    std::string message_m;
  };
  std::vector<PhaseDiagnostic> phase_diagnostics_m;
};

struct TunaSimulationCase {
  int simulation_id_m = 0;
  TunaAssessmentRunSummary truth_m;
  TunaAssessmentRunSummary estimated_m;
  double depletion_bias_m = std::numeric_limits<double>::quiet_NaN();
};

struct TunaSimulationResult {
  std::vector<TunaSimulationCase> cases_m;
  double mean_depletion_bias_m = std::numeric_limits<double>::quiet_NaN();
  double median_depletion_bias_m = std::numeric_limits<double>::quiet_NaN();
  double p10_depletion_bias_m = std::numeric_limits<double>::quiet_NaN();
  double p90_depletion_bias_m = std::numeric_limits<double>::quiet_NaN();
  double low_depletion_threshold_m = 0.05;
  double low_depletion_rate_m = std::numeric_limits<double>::quiet_NaN();
  int low_depletion_count_m = 0;
  int n_finite_bias_m = 0;
};

inline double empirical_quantile(const std::vector<double> &sorted_values,
                                 double probability) {
  if (sorted_values.empty()) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  const double p = std::max(0.0, std::min(1.0, probability));
  const double idx = p * static_cast<double>(sorted_values.size() - 1);
  const size_t lo = static_cast<size_t>(std::floor(idx));
  const size_t hi = static_cast<size_t>(std::ceil(idx));

  if (lo == hi) {
    return sorted_values[lo];
  }

  const double w = idx - static_cast<double>(lo);
  return sorted_values[lo] + w * (sorted_values[hi] - sorted_values[lo]);
}

inline double report_value_or_nan(const ModelReportContext &ctx,
                                  const std::string &name) {
  for (const auto &rv : ctx.reports().values()) {
    if (rv.name_m == name) {
      return rv.value_m;
    }
  }
  return std::numeric_limits<double>::quiet_NaN();
}

inline TunaAssessmentRunSummary
evaluate_at_parameters(const TunaSpatialAssessmentData &data,
                       const TunaAssessmentControls &controls,
                       const std::vector<double> &parameters,
                       const std::string &label = "base") {
  AdvancedSpatialTunaAssessmentModel model(data, controls);
  const ParameterSet pset = model.parameter_set();

  if (parameters.size() != pset.size()) {
    throw std::invalid_argument(
        "evaluate_at_parameters: parameter length mismatch");
  }

  ModelReportContext ctx;
  const double nll = model.evaluate(parameters, ctx);

  TunaAssessmentRunSummary out;
  out.label_m = label;
  out.n_years_m = data.n_years_m;
  out.nll_m = nll;
  out.r0_m = report_value_or_nan(ctx, "r0");
  out.steepness_m = report_value_or_nan(ctx, "steepness");
  out.ssb0_m = report_value_or_nan(ctx, "ssb0");
  out.ssb_terminal_m = report_value_or_nan(ctx, "ssb_terminal");
  out.depletion_terminal_m = report_value_or_nan(ctx, "depletion_terminal");

  out.ssb_by_year_m.assign(static_cast<size_t>(data.n_years_m),
                           std::numeric_limits<double>::quiet_NaN());
  for (int y = 0; y < data.n_years_m; ++y) {
    out.ssb_by_year_m[static_cast<size_t>(y)] =
        report_value_or_nan(ctx, "ssb_year_" + std::to_string(y + 1));
  }

  return out;
}

inline std::string diagnostic_row_key(int y, int s, int f, int r) {
  return "y" + std::to_string(y + 1) + "_s" + std::to_string(s + 1) + "_f" +
         std::to_string(f + 1) + "_r" + std::to_string(r + 1);
}

inline double sample_sd(const std::vector<double> &values, double mean) {
  if (values.size() < 2) {
    return std::numeric_limits<double>::quiet_NaN();
  }
  double ss = 0.0;
  for (double value : values) {
    const double d = value - mean;
    ss += d * d;
  }
  return std::sqrt(ss / static_cast<double>(values.size() - 1));
}

inline TunaDiagnosticBundle
evaluate_diagnostic_bundle(const TunaSpatialAssessmentData &data,
                           const TunaAssessmentControls &controls,
                           const std::vector<double> &parameters) {
  TunaAssessmentControls diagnostic_controls = controls;
  diagnostic_controls.phase_m = TunaAssessmentPhase::Full;
  diagnostic_controls.report_observation_predictions_m = true;
  AdvancedSpatialTunaAssessmentModel model(data, diagnostic_controls);
  if (parameters.size() != model.parameter_set().size()) {
    throw std::invalid_argument(
        "evaluate_diagnostic_bundle: parameter length mismatch");
  }

  ModelReportContext ctx;
  model.evaluate(parameters, ctx);
  TunaDiagnosticBundle out;

  const std::vector<std::string> likelihood_names = {
      "prior",   "availability_penalty", "movement_penalty",
      "index",   "composition",          "retained",
      "discard", "recruitment"};
  double absolute_total = 0.0;
  for (const std::string &name : likelihood_names) {
    TunaLikelihoodComponent row;
    row.component_m = name;
    row.nll_m = report_value_or_nan(ctx, "nll_" + name);
    if (std::isfinite(row.nll_m))
      absolute_total += std::abs(row.nll_m);
    out.likelihood_components_m.push_back(row);
  }
  for (auto &row : out.likelihood_components_m) {
    if (row.component_m == "prior")
      row.observation_units_m = static_cast<int>(parameters.size());
    else if (row.component_m == "recruitment")
      row.observation_units_m = data.n_years_m;
    else if (row.component_m == "availability_penalty" ||
             row.component_m == "movement_penalty")
      row.observation_units_m = 0;
    else {
      for (int f = 0; f < data.n_fleets_m; ++f)
        for (int y = 0; y < data.n_years_m; ++y)
          for (int s = 0; s < data.n_seasons_m; ++s)
            for (int r = 0; r < data.n_regions_m; ++r) {
              const size_t idx =
                  data.fleet_year_season_region_index(f, y, s, r);
              if (row.component_m == "index" &&
                  diagnostic_controls.use_index_likelihood_m &&
                  data.observed_index_m[idx] > 0.0)
                ++row.observation_units_m;
              else if (row.component_m == "retained" &&
                       diagnostic_controls.use_retained_biomass_likelihood_m &&
                       data.observed_retained_biomass_m[idx] > 0.0)
                ++row.observation_units_m;
              else if (row.component_m == "discard" &&
                       diagnostic_controls.use_discard_biomass_likelihood_m &&
                       data.observed_discard_biomass_m[idx] > 0.0)
                ++row.observation_units_m;
              else if (row.component_m == "composition" &&
                       diagnostic_controls.use_catch_composition_likelihood_m) {
                int total = 0;
                for (int a = 0; a < data.n_ages_m; ++a)
                  total += data.observed_catch_numbers_m
                               [data.fleet_year_season_region_age_index(f, y, s,
                                                                        r, a)];
                if (total > 0)
                  ++row.observation_units_m;
              }
            }
    }
    if (row.observation_units_m > 0)
      row.nll_per_unit_m =
          row.nll_m / static_cast<double>(row.observation_units_m);
  }
  if (absolute_total > 0.0) {
    for (auto &row : out.likelihood_components_m)
      row.share_absolute_m = std::abs(row.nll_m) / absolute_total;
  }

  const ParameterSet parameter_set = model.parameter_set();
  for (size_t i = 0; i < parameters.size(); ++i) {
    TunaParameterDiagnostic row;
    const ParameterInfo &info = parameter_set.values()[i];
    row.name_m = info.name_m;
    row.value_m = parameters[i];
    row.initial_m = info.initial_value_m;
    row.displacement_m = row.value_m - row.initial_m;
    row.random_effect_m = info.is_random_m;
    row.extreme_m = !std::isfinite(row.value_m) ||
                    std::abs(row.value_m) > 20.0 ||
                    std::abs(row.displacement_m) > 8.0;
    out.parameters_m.push_back(row);
  }

  for (int f = 0; f < data.n_fleets_m; ++f) {
    for (int s = 0; s < data.n_seasons_m; ++s) {
      for (int r = 0; r < data.n_regions_m; ++r) {
        TunaCatchabilityAvailabilityDiagnostic row;
        row.fleet_m = f + 1;
        row.season_m = s + 1;
        row.region_m = r + 1;
        row.log_fishing_q_m =
            report_value_or_nan(ctx, "diag_log_q_f" + std::to_string(f + 1));
        row.log_index_q_m = report_value_or_nan(ctx, "diag_log_index_q_f" +
                                                         std::to_string(f + 1));
        row.log_availability_m = report_value_or_nan(
            ctx, "diag_log_availability_f" + std::to_string(f + 1) + "_s" +
                     std::to_string(s + 1) + "_r" + std::to_string(r + 1));
        row.log_fishing_q_availability_m = report_value_or_nan(
            ctx, "diag_log_q_availability_f" + std::to_string(f + 1) + "_s" +
                     std::to_string(s + 1) + "_r" + std::to_string(r + 1));
        out.catchability_availability_m.push_back(row);
      }
    }
  }

  const auto add_lognormal = [&](const std::string &component, double observed,
                                 int y, int s, int f, int r) {
    if (!(observed > 0.0)) {
      return;
    }
    const std::string prefix = "diag_" + component + "_";
    const std::string key = diagnostic_row_key(y, s, f, r);
    const double predicted = report_value_or_nan(ctx, prefix + "pred_" + key);
    const double sigma = report_value_or_nan(ctx, prefix + "sigma_" + key);
    TunaObservationDiagnostic row;
    row.component_m = component;
    row.year_m = y + 1;
    row.season_m = s + 1;
    row.fleet_m = f + 1;
    row.region_m = r + 1;
    row.observed_m = observed;
    row.predicted_m = predicted;
    row.standard_deviation_m = sigma;
    if (predicted > 0.0 && sigma > 0.0) {
      row.standardized_residual_m =
          (std::log(observed) - std::log(predicted)) / sigma;
    }
    out.observations_m.push_back(row);
  };

  for (int f = 0; f < data.n_fleets_m; ++f) {
    for (int y = 0; y < data.n_years_m; ++y) {
      for (int s = 0; s < data.n_seasons_m; ++s) {
        for (int r = 0; r < data.n_regions_m; ++r) {
          const size_t idx = data.fleet_year_season_region_index(f, y, s, r);
          if (diagnostic_controls.use_index_likelihood_m)
            add_lognormal("index", data.observed_index_m[idx], y, s, f, r);
          if (diagnostic_controls.use_retained_biomass_likelihood_m)
            add_lognormal("retained", data.observed_retained_biomass_m[idx], y,
                          s, f, r);
          if (diagnostic_controls.use_discard_biomass_likelihood_m)
            add_lognormal("discard", data.observed_discard_biomass_m[idx], y, s,
                          f, r);

          int total = 0;
          for (int a = 0; a < data.n_ages_m; ++a)
            total +=
                data.observed_catch_numbers_m
                    [data.fleet_year_season_region_age_index(f, y, s, r, a)];
          if (!diagnostic_controls.use_catch_composition_likelihood_m ||
              total <= 0)
            continue;

          const std::string key = diagnostic_row_key(y, s, f, r);
          const double theta =
              report_value_or_nan(ctx, "diag_comp_theta_" + key);
          const double variance_multiplier =
              1.0 + static_cast<double>(total - 1) * theta / (1.0 + theta);
          for (int a = 0; a < data.n_ages_m; ++a) {
            TunaObservationDiagnostic row;
            row.component_m = "composition";
            row.year_m = y + 1;
            row.season_m = s + 1;
            row.fleet_m = f + 1;
            row.region_m = r + 1;
            row.age_m = a + 1;
            const int observed_count =
                data.observed_catch_numbers_m
                    [data.fleet_year_season_region_age_index(f, y, s, r, a)];
            const double p = report_value_or_nan(
                ctx, "diag_comp_pred_" + key + "_a" + std::to_string(a + 1));
            row.observed_m = static_cast<double>(observed_count) /
                             static_cast<double>(total);
            row.predicted_m = p;
            if (p >= 0.0 && p <= 1.0 && theta >= 0.0) {
              const double count_variance = static_cast<double>(total) * p *
                                            (1.0 - p) * variance_multiplier;
              row.standard_deviation_m =
                  std::sqrt(std::max(0.0, count_variance)) /
                  static_cast<double>(total);
              if (row.standard_deviation_m > 0.0) {
                row.standardized_residual_m =
                    (row.observed_m - p) / row.standard_deviation_m;
              }
            }
            out.observations_m.push_back(row);
          }
        }
      }
    }
  }

  const std::vector<std::string> components = {"index", "retained", "discard",
                                               "composition"};
  for (const std::string &component : components) {
    for (int fleet = 1; fleet <= data.n_fleets_m; ++fleet) {
      std::vector<double> residuals;
      double squared_error = 0.0;
      int prediction_count = 0;
      for (const auto &row : out.observations_m) {
        if (row.component_m != component || row.fleet_m != fleet)
          continue;
        if (std::isfinite(row.standardized_residual_m))
          residuals.push_back(row.standardized_residual_m);
        if (std::isfinite(row.observed_m) && std::isfinite(row.predicted_m)) {
          const double error = row.observed_m - row.predicted_m;
          squared_error += error * error;
          ++prediction_count;
        }
      }
      if (residuals.empty() && prediction_count == 0)
        continue;
      TunaResidualSummary summary;
      summary.component_m = component;
      summary.fleet_m = fleet;
      summary.n_m = static_cast<int>(residuals.size());
      if (!residuals.empty()) {
        double sum = 0.0;
        for (double residual : residuals)
          sum += residual;
        summary.mean_residual_m = sum / static_cast<double>(residuals.size());
        summary.sdnr_m = sample_sd(residuals, summary.mean_residual_m);
      }
      if (prediction_count > 0)
        summary.rmse_m = std::sqrt(squared_error / prediction_count);
      out.residual_summaries_m.push_back(summary);
    }
  }

  const std::vector<std::string> strata = {"year", "season", "region"};
  for (const std::string &component : components) {
    for (int fleet = 1; fleet <= data.n_fleets_m; ++fleet) {
      for (const std::string &stratum : strata) {
        const int levels = stratum == "year"     ? data.n_years_m
                           : stratum == "season" ? data.n_seasons_m
                                                 : data.n_regions_m;
        for (int level = 1; level <= levels; ++level) {
          std::vector<double> residuals;
          double squared_error = 0.0;
          int prediction_count = 0;
          for (const auto &row : out.observations_m) {
            const int row_level = stratum == "year"     ? row.year_m
                                  : stratum == "season" ? row.season_m
                                                        : row.region_m;
            if (row.component_m != component || row.fleet_m != fleet ||
                row_level != level)
              continue;
            if (std::isfinite(row.standardized_residual_m))
              residuals.push_back(row.standardized_residual_m);
            if (std::isfinite(row.observed_m) &&
                std::isfinite(row.predicted_m)) {
              const double error = row.observed_m - row.predicted_m;
              squared_error += error * error;
              ++prediction_count;
            }
          }
          if (residuals.empty() && prediction_count == 0)
            continue;
          TunaStratifiedResidualSummary summary;
          summary.component_m = component;
          summary.fleet_m = fleet;
          summary.stratum_m = stratum;
          summary.level_m = level;
          summary.n_m = static_cast<int>(residuals.size());
          if (!residuals.empty()) {
            double sum = 0.0;
            for (double residual : residuals)
              sum += residual;
            summary.mean_residual_m = sum / residuals.size();
            summary.sdnr_m = sample_sd(residuals, summary.mean_residual_m);
          }
          if (prediction_count > 0)
            summary.rmse_m = std::sqrt(squared_error / prediction_count);
          out.stratified_residuals_m.push_back(summary);
        }
      }
    }
  }
  return out;
}

inline TunaAssessmentRunSummary
evaluate_initial_run(const TunaSpatialAssessmentData &data,
                     const TunaAssessmentControls &controls,
                     const std::string &label = "base") {
  AdvancedSpatialTunaAssessmentModel model(data, controls);
  return evaluate_at_parameters(data, controls,
                                model.parameter_set().initials(), label);
}

inline bool has_prefix(const std::string &value, const std::string &prefix) {
  return value.size() >= prefix.size() &&
         value.compare(0, prefix.size(), prefix) == 0;
}

inline std::vector<bool>
build_phase_free_mask(const std::vector<std::string> &fixed_names,
                      TunaAssessmentPhase phase) {
  std::vector<bool> free_mask(fixed_names.size(), false);

  for (size_t i = 0; i < fixed_names.size(); ++i) {
    const std::string &name = fixed_names[i];

    const bool is_recruit_core = name == "log_r0" ||
                                 name == "logit_steepness" ||
                                 name == "log_sigma_recruit";

    const bool is_catchability_block =
        has_prefix(name, "log_q_fleet_") ||
        has_prefix(name, "log_index_q_fleet_") ||
        has_prefix(name, "sel50_raw_fleet_") ||
        has_prefix(name, "log_sel_slope_fleet_") ||
        has_prefix(name, "retention50_raw_fleet_") ||
        has_prefix(name, "log_retention_slope_fleet_") ||
        has_prefix(name, "log_sigma_index_fleet_") ||
        has_prefix(name, "log_theta_comp_fleet_") ||
        has_prefix(name, "log_sigma_retained_bio_fleet_") ||
        has_prefix(name, "log_sigma_discard_bio_fleet_");

    const bool is_movement_block = has_prefix(name, "move_logit_season_");
    const bool is_availability_block =
        has_prefix(name, "log_availability_scale_fleet_");

    if (phase == TunaAssessmentPhase::InitializeRecruitment) {
      free_mask[i] = is_recruit_core;
    } else if (phase == TunaAssessmentPhase::InitializeCatchability) {
      free_mask[i] = is_recruit_core || is_catchability_block;
    } else if (phase == TunaAssessmentPhase::InitializeMovement) {
      free_mask[i] = is_recruit_core || is_catchability_block ||
                     is_movement_block || is_availability_block;
    } else {
      free_mask[i] = true;
    }
  }

  return free_mask;
}

inline int parse_positive_int_after_key(const std::string &text,
                                        const std::string &key) {
  const size_t pos = text.find(key);
  if (pos == std::string::npos) {
    return -1;
  }

  size_t i = pos + key.size();
  if (i >= text.size() || text[i] < '0' || text[i] > '9') {
    return -1;
  }

  int value = 0;
  while (i < text.size() && text[i] >= '0' && text[i] <= '9') {
    value = value * 10 + static_cast<int>(text[i] - '0');
    ++i;
  }

  return value;
}

inline std::vector<bool>
build_partial_movement_free_mask(const std::vector<std::string> &fixed_names) {
  std::vector<bool> free_mask(fixed_names.size(), false);

  for (size_t i = 0; i < fixed_names.size(); ++i) {
    const std::string &name = fixed_names[i];

    const bool is_recruit_core = name == "log_r0" ||
                                 name == "logit_steepness" ||
                                 name == "log_sigma_recruit";

    const bool is_catchability_block =
        has_prefix(name, "log_q_fleet_") ||
        has_prefix(name, "log_index_q_fleet_") ||
        has_prefix(name, "sel50_raw_fleet_") ||
        has_prefix(name, "log_sel_slope_fleet_") ||
        has_prefix(name, "retention50_raw_fleet_") ||
        has_prefix(name, "log_retention_slope_fleet_") ||
        has_prefix(name, "log_sigma_index_fleet_") ||
        has_prefix(name, "log_theta_comp_fleet_") ||
        has_prefix(name, "log_sigma_retained_bio_fleet_") ||
        has_prefix(name, "log_sigma_discard_bio_fleet_");

    bool is_diag_movement = false;
    if (has_prefix(name, "move_logit_season_")) {
      const int from_idx = parse_positive_int_after_key(name, "_from_");
      const int to_idx = parse_positive_int_after_key(name, "_to_");
      if (from_idx > 0 && to_idx > 0) {
        // For an n-1 logit parameterization, to==from encodes staying.
        is_diag_movement = (from_idx == to_idx);
      }
    }

    free_mask[i] = is_recruit_core || is_catchability_block || is_diag_movement;
  }

  return free_mask;
}

inline std::vector<size_t> build_movement_availability_full_indices(
    const std::vector<std::string> &fixed_names,
    const ParameterPartition &partition) {
  std::vector<size_t> indices;
  indices.reserve(fixed_names.size());

  for (size_t i = 0; i < fixed_names.size(); ++i) {
    const std::string &name = fixed_names[i];
    if (has_prefix(name, "move_logit_season_") ||
        has_prefix(name, "log_availability_scale_fleet_")) {
      indices.push_back(partition.fixed_indices_m[i]);
    }
  }

  return indices;
}

struct FixedTrustRegionSpec {
  std::vector<size_t> full_indices_m;
  std::vector<double> centers_m;
  std::vector<double> half_widths_m;
};

inline FixedTrustRegionSpec build_movement_entry_trust_region_spec(
    const std::vector<std::string> &fixed_names,
    const ParameterPartition &partition, const std::vector<double> &theta_fixed,
    double movement_half_width = 0.35, double availability_half_width = 0.25) {
  if (theta_fixed.size() != fixed_names.size()) {
    throw std::invalid_argument("build_movement_entry_trust_region_spec: "
                                "theta/fixed-name length mismatch");
  }

  FixedTrustRegionSpec spec;
  spec.full_indices_m.reserve(fixed_names.size());
  spec.centers_m.reserve(fixed_names.size());
  spec.half_widths_m.reserve(fixed_names.size());

  for (size_t i = 0; i < fixed_names.size(); ++i) {
    const std::string &name = fixed_names[i];
    if (has_prefix(name, "move_logit_season_")) {
      spec.full_indices_m.push_back(partition.fixed_indices_m[i]);
      spec.centers_m.push_back(theta_fixed[i]);
      spec.half_widths_m.push_back(movement_half_width);
    } else if (has_prefix(name, "log_availability_scale_fleet_")) {
      spec.full_indices_m.push_back(partition.fixed_indices_m[i]);
      spec.centers_m.push_back(theta_fixed[i]);
      spec.half_widths_m.push_back(availability_half_width);
    }
  }

  return spec;
}

inline std::vector<double>
build_start_jitter_scales(const std::vector<std::string> &fixed_names) {
  std::vector<double> scales(fixed_names.size(), 0.20);

  for (size_t i = 0; i < fixed_names.size(); ++i) {
    const std::string &name = fixed_names[i];

    if (name == "log_r0" || name == "logit_steepness" ||
        name == "log_sigma_recruit") {
      scales[i] = 0.30;
    } else if (has_prefix(name, "log_q_fleet_") ||
               has_prefix(name, "log_index_q_fleet_") ||
               has_prefix(name, "sel50_raw_fleet_") ||
               has_prefix(name, "log_sel_slope_fleet_") ||
               has_prefix(name, "retention50_raw_fleet_") ||
               has_prefix(name, "log_retention_slope_fleet_")) {
      scales[i] = 0.18;
    } else if (has_prefix(name, "move_logit_season_") ||
               has_prefix(name, "log_availability_scale_fleet_")) {
      scales[i] = 0.10;
    } else if (has_prefix(name, "log_sigma_") ||
               has_prefix(name, "log_theta_comp_fleet_")) {
      scales[i] = 0.12;
    }
  }

  return scales;
}

inline const char *phase_name(TunaAssessmentPhase phase) {
  switch (phase) {
  case TunaAssessmentPhase::InitializeRecruitment:
    return "InitializeRecruitment";
  case TunaAssessmentPhase::InitializeCatchability:
    return "InitializeCatchability";
  case TunaAssessmentPhase::InitializeMovement:
    return "InitializeMovement";
  case TunaAssessmentPhase::Full:
    return "Full";
  }

  return "Unknown";
}

inline bool
is_initial_exact_gradient_failure(const LaplaceExactLBFGSResult &fit) {
  return fit.iterations_m == 0 && !fit.logdet_ok_m &&
         fit.message_m.find(
             "initial exact Laplace gradient evaluation failed") !=
             std::string::npos;
}

inline TunaAssessmentControls
tuned_controls_for_phase(const TunaAssessmentControls &base,
                         TunaAssessmentPhase phase) {
  TunaAssessmentControls controls = base;
  controls.phase_m = phase;

  if (phase == TunaAssessmentPhase::InitializeRecruitment) {
    controls.use_priors_m = true;
    controls.sd_prior_log_q_m = std::min(controls.sd_prior_log_q_m, 0.25);
    controls.sd_prior_log_sel50_m =
        std::min(controls.sd_prior_log_sel50_m, 0.30);
    controls.sd_prior_log_sel_slope_m =
        std::min(controls.sd_prior_log_sel_slope_m, 0.30);
    controls.sd_prior_log_availability_scale_m =
        std::min(controls.sd_prior_log_availability_scale_m, 0.20);
    controls.sd_prior_move_logit_m =
        std::min(controls.sd_prior_move_logit_m, 0.20);
    controls.sd_prior_log_sigma_m =
        std::min(controls.sd_prior_log_sigma_m, 0.60);
    controls.movement_smoothing_weight_m =
        std::max(controls.movement_smoothing_weight_m, 35.0);
    controls.availability_smoothing_weight_m =
        std::max(controls.availability_smoothing_weight_m, 35.0);
  } else if (phase == TunaAssessmentPhase::InitializeCatchability) {
    controls.use_priors_m = true;
    controls.sd_prior_log_q_m = std::min(controls.sd_prior_log_q_m, 0.55);
    controls.sd_prior_log_sel50_m =
        std::min(controls.sd_prior_log_sel50_m, 0.65);
    controls.sd_prior_log_sel_slope_m =
        std::min(controls.sd_prior_log_sel_slope_m, 0.65);
    controls.sd_prior_log_availability_scale_m =
        std::min(controls.sd_prior_log_availability_scale_m, 0.45);
    controls.sd_prior_move_logit_m =
        std::min(controls.sd_prior_move_logit_m, 0.45);
    controls.movement_smoothing_weight_m =
        std::max(controls.movement_smoothing_weight_m, 20.0);
    controls.availability_smoothing_weight_m =
        std::max(controls.availability_smoothing_weight_m, 20.0);
  } else if (phase == TunaAssessmentPhase::InitializeMovement) {
    controls.use_priors_m = true;
    controls.sd_prior_log_availability_scale_m =
        std::min(controls.sd_prior_log_availability_scale_m, 0.80);
    controls.sd_prior_move_logit_m =
        std::min(controls.sd_prior_move_logit_m, 0.80);
    controls.movement_smoothing_weight_m =
        std::max(controls.movement_smoothing_weight_m, 10.0);
    controls.availability_smoothing_weight_m =
        std::max(controls.availability_smoothing_weight_m, 10.0);
  }

  return controls;
}

class PhaseLockedSpatialModel : public QuadraModel<PhaseLockedSpatialModel> {
public:
  PhaseLockedSpatialModel(TunaSpatialAssessmentData data,
                          TunaAssessmentControls controls,
                          std::vector<size_t> locked_full_indices,
                          std::vector<double> locked_full_values,
                          std::vector<size_t> continuation_full_indices = {},
                          double continuation_factor = 1.0,
                          std::vector<size_t> trust_full_indices = {},
                          std::vector<double> trust_centers = {},
                          std::vector<double> trust_half_widths = {})
      : base_model_m(std::move(data), controls),
        locked_full_indices_m(std::move(locked_full_indices)),
        locked_full_values_m(std::move(locked_full_values)),
        continuation_full_indices_m(std::move(continuation_full_indices)),
        continuation_factor_m(continuation_factor),
        trust_full_indices_m(std::move(trust_full_indices)),
        trust_centers_m(std::move(trust_centers)),
        trust_half_widths_m(std::move(trust_half_widths)) {
    if (locked_full_indices_m.size() != locked_full_values_m.size()) {
      throw std::invalid_argument(
          "PhaseLockedSpatialModel: lock index/value length mismatch");
    }

    if (!(continuation_factor_m > 0.0 && continuation_factor_m <= 1.0)) {
      throw std::invalid_argument(
          "PhaseLockedSpatialModel: continuation factor must be in (0, 1]");
    }

    if (trust_full_indices_m.size() != trust_centers_m.size() ||
        trust_full_indices_m.size() != trust_half_widths_m.size()) {
      throw std::invalid_argument(
          "PhaseLockedSpatialModel: trust-region vectors length mismatch");
    }
  }

  std::vector<std::string> parameter_names_impl() const {
    return base_model_m.parameter_names_impl();
  }

  template <typename Type>
  Type evaluate_impl(const std::vector<Type> &parameters,
                     ModelReportContext &ctx) const {
    std::vector<Type> projected = parameters;
    for (size_t i = 0; i < locked_full_indices_m.size(); ++i) {
      const size_t idx = locked_full_indices_m[i];
      projected[idx] = Type(locked_full_values_m[i]);
    }

    if (continuation_factor_m < 1.0) {
      const Type cf = Type(continuation_factor_m);
      for (size_t idx : continuation_full_indices_m) {
        projected[idx] = cf * projected[idx];
      }
    }

    for (size_t i = 0; i < trust_full_indices_m.size(); ++i) {
      const size_t idx = trust_full_indices_m[i];
      const Type center = Type(trust_centers_m[i]);
      const Type half_width = Type(trust_half_widths_m[i]);
      const Type lo = center - half_width;
      const Type hi = center + half_width;

      if constexpr (std::is_same_v<Type, had::ThirdOrderScalar>) {
        if (projected[idx].val < lo.val)
          projected[idx] = lo;
        else if (projected[idx].val > hi.val)
          projected[idx] = hi;
      } else {
        if (projected[idx] < lo)
          projected[idx] = lo;
        else if (projected[idx] > hi)
          projected[idx] = hi;
      }
    }

    return base_model_m.template evaluate<Type>(projected, ctx);
  }

private:
  AdvancedSpatialTunaAssessmentModel base_model_m;
  std::vector<size_t> locked_full_indices_m;
  std::vector<double> locked_full_values_m;
  std::vector<size_t> continuation_full_indices_m;
  double continuation_factor_m = 1.0;
  std::vector<size_t> trust_full_indices_m;
  std::vector<double> trust_centers_m;
  std::vector<double> trust_half_widths_m;
};

inline LaplaceExactLBFGSOptions
make_exact_lbfgs_options(const TunaFitOptions &options) {
  LaplaceExactLBFGSOptions out;
  out.max_iterations_m = options.max_iterations_per_phase_m;
  out.memory_m = options.lbfgs_memory_m;
  out.max_evaluations_m = options.max_evaluations_per_phase_m;
  out.gradient_tolerance_m = options.gradient_tolerance_m;
  out.initial_step_scale_m = options.initial_step_m;
  out.min_step_scale_m = options.min_step_m;
  out.use_backtracking_m = true;
  out.warm_start_random_m = true;
  out.print_every_m = options.lbfgs_print_every_m;
  return out;
}

template <class Model>
inline LaplaceExactLBFGSResult optimize_full_phase_public_exact_laplace(
    Model &model, const std::vector<double> &theta,
    const std::vector<double> &random, const ParameterPartition &partition,
    const TunaFitOptions &options,
    const std::vector<std::string> &fixed_names) {
  LaplaceObjectiveOptions objective_options;
  objective_options.include_constant_m = true;
  objective_options.compute_mixed_derivatives_m = true;
  objective_options.newton_m.max_iterations_m =
      std::max(300, options.max_iterations_per_phase_m + 50);
  objective_options.newton_m.gradient_tolerance_m = 1e-5;
  objective_options.newton_m.step_tolerance_m = 1e-10;
  objective_options.newton_m.initial_step_scale_m = 1.0;
  objective_options.newton_m.min_step_scale_m = 1e-10;
  objective_options.newton_m.sufficient_decrease_m = 1e-6;
  objective_options.newton_m.hessian_drop_tol_m = 0.0;
  objective_options.newton_m.use_backtracking_m = true;

  laplace::ExactLaplaceGradientEngineOptions engine_options;
  engine_options.discover_active_directions = false;
  engine_options.stream_dense_hdot_trace = true;
  engine_options.hdot_workers = options.hdot_workers_m;

  stats::ExactLaplaceEvaluator<Model> evaluator(
      model, theta, random, partition, objective_options, engine_options);

  stats::LaplaceOptimizerOptions public_options;
  public_options.max_iterations = options.max_iterations_per_phase_m;
  public_options.memory = std::max(20, options.lbfgs_memory_m);
  public_options.max_evaluations = options.max_evaluations_per_phase_m;
  public_options.gradient_tolerance = options.gradient_tolerance_m;
  public_options.step_tolerance = 1e-10;
  public_options.initial_step_scale = 1.0;
  public_options.maximum_direction_norm = 1.0;
  public_options.minimum_step_scale = std::min(options.min_step_m, 1e-10);
  public_options.sufficient_decrease = 1e-4;
  public_options.print_every = options.lbfgs_print_every_m;
  public_options.gradient_table_every = 10;
  public_options.gradient_names = fixed_names;

  const stats::LaplaceOptimizerResult public_fit =
      stats::optimize_laplace(evaluator, theta, public_options);

  LaplaceExactLBFGSResult out;
  out.theta_initial_m = theta;
  out.theta_hat_m = public_fit.fixed;
  out.u_hat_m = public_fit.random_mode;
  out.full_hat_m = public_fit.full;
  out.gradient_fixed_m = public_fit.gradient;
  out.laplace_objective_m = public_fit.objective;
  out.gradient_norm_m = public_fit.gradient_norm;
  out.step_norm_m = public_fit.step_norm;
  out.iterations_m = public_fit.iterations;
  out.converged_m = public_fit.converged;
  out.logdet_ok_m = public_fit.evaluation.objective.logdet_ok_m;
  out.message_m = "public_exact_laplace: " + public_fit.message;
  const bool near_stationary_line_search_stop =
      !out.converged_m && out.logdet_ok_m &&
      std::isfinite(out.gradient_norm_m) &&
      out.gradient_norm_m <=
          options.line_search_acceptance_gradient_tolerance_m &&
      public_fit.message.find("line search") != std::string::npos;
  if (near_stationary_line_search_stop) {
    out.converged_m = true;
    out.message_m +=
        " Accepted as near-stationary: verified exact Laplace gradient "
        "is below the configured line-search acceptance tolerance.";
  }
  for (const auto &item : public_fit.history) {
    LaplaceLBFGSIteration converted;
    converted.iteration_m = item.iteration;
    converted.objective_m = item.objective;
    converted.gradient_norm_m = item.gradient_norm;
    converted.step_scale_m = item.step_scale;
    converted.step_norm_m = item.step_norm;
    out.history_m.push_back(converted);
  }
  return out;
}

inline JointOnlyExactLBFGSOptions
make_joint_only_lbfgs_options(const TunaFitOptions &options) {
  JointOnlyExactLBFGSOptions out;
  out.max_iterations_m = options.max_iterations_per_phase_m;
  out.memory_m = options.lbfgs_memory_m;
  out.max_evaluations_m = options.max_evaluations_per_phase_m;
  out.print_every_m = options.lbfgs_print_every_m;
  out.gradient_tolerance_m = options.gradient_tolerance_m;
  out.initial_step_scale_m = options.initial_step_m;
  out.min_step_scale_m = options.min_step_m;
  out.use_backtracking_m = true;
  out.warm_start_random_m = true;
  return out;
}

inline JointOnlyExactLBFGSOptions
make_robust_joint_only_lbfgs_options(const TunaFitOptions &options) {
  JointOnlyExactLBFGSOptions out = make_joint_only_lbfgs_options(options);

  out.max_iterations_m =
      std::max(20, std::min(options.max_iterations_per_phase_m + 10, 80));
  out.gradient_tolerance_m = std::max(options.gradient_tolerance_m, 1e-4);
  out.initial_step_scale_m = std::min(out.initial_step_scale_m, 0.02);
  out.min_step_scale_m = std::min(out.min_step_scale_m, 1e-10);
  out.warm_start_random_m = false;

  out.gradient_m.newton_m.max_iterations_m =
      std::max(40, options.max_iterations_per_phase_m + 10);
  out.gradient_m.newton_m.gradient_tolerance_m = 1e-6;
  out.gradient_m.newton_m.step_tolerance_m = 1e-10;
  out.gradient_m.newton_m.initial_step_scale_m = 0.25;
  out.gradient_m.newton_m.min_step_scale_m = 1e-10;
  out.gradient_m.newton_m.sufficient_decrease_m = 1e-6;
  out.gradient_m.newton_m.hessian_drop_tol_m = 1e-10;
  out.gradient_m.newton_m.use_backtracking_m = true;

  return out;
}

inline RandomEffectNewtonOptions
make_movement_warmup_newton_options(const TunaFitOptions &options) {
  RandomEffectNewtonOptions out;
  out.max_iterations_m = std::max(40, options.max_iterations_per_phase_m + 15);
  out.gradient_tolerance_m = 1e-6;
  out.step_tolerance_m = 1e-10;
  out.initial_step_scale_m = 0.20;
  out.min_step_scale_m = 1e-10;
  out.sufficient_decrease_m = 1e-6;
  out.hessian_drop_tol_m = 1e-10;
  out.use_backtracking_m = true;
  return out;
}

inline TunaFitResult
fit_spatial_assessment(const TunaSpatialAssessmentData &data,
                       const TunaAssessmentControls &base_controls,
                       const TunaFitOptions &options) {
  AdvancedSpatialTunaAssessmentModel prototype(data, base_controls);
  const ParameterSet pset = prototype.parameter_set();
  const ParameterPartition partition = partition_parameters(pset);
  const std::vector<std::string> fixed_names =
      parameter_names_by_indices(pset, partition.fixed_indices_m);
  const std::vector<size_t> movement_availability_full_indices =
      build_movement_availability_full_indices(fixed_names, partition);
  const std::vector<double> jitter_scales =
      build_start_jitter_scales(fixed_names);
  const std::vector<double> full_initial = pset.initials();
  const PartitionedVector<double> split0 =
      split_parameters(full_initial, partition);

  if (options.multistart_m <= 0) {
    throw std::invalid_argument(
        "fit_spatial_assessment: multistart must be positive");
  }
  if (options.fixed_parameter_names_m.size() !=
      options.fixed_parameter_values_m.size()) {
    throw std::invalid_argument(
        "fit_spatial_assessment: fixed parameter names/values mismatch");
  }
  if (!std::isfinite(options.line_search_acceptance_gradient_tolerance_m) ||
      options.line_search_acceptance_gradient_tolerance_m < 0.0) {
    throw std::invalid_argument(
        "fit_spatial_assessment: line-search acceptance gradient "
        "tolerance must be finite and non-negative");
  }

  TunaFitResult best_fit;
  double best_obj = std::numeric_limits<double>::infinity();

  std::mt19937_64 rng(options.seed_m);
  std::normal_distribution<double> jitter(0.0, options.jitter_sd_m);

  for (int start_idx = 0; start_idx < options.multistart_m; ++start_idx) {
    if (options.lbfgs_print_every_m > 0) {
      std::cout << "tuna_fit start=" << (start_idx + 1) << '/'
                << options.multistart_m << " begin\n"
                << std::flush;
    }
    std::vector<double> theta = split0.fixed_m;
    std::vector<double> random = split0.random_m;
    if (start_idx > 0) {
      for (size_t j = 0; j < theta.size(); ++j) {
        theta[j] += jitter(rng) * jitter_scales[j];
      }
    }
    for (size_t lock = 0; lock < options.fixed_parameter_names_m.size();
         ++lock) {
      const auto it = std::find(fixed_names.begin(), fixed_names.end(),
                                options.fixed_parameter_names_m[lock]);
      if (it == fixed_names.end()) {
        throw std::invalid_argument(
            "fit_spatial_assessment: unknown fixed parameter " +
            options.fixed_parameter_names_m[lock]);
      }
      theta[static_cast<size_t>(std::distance(fixed_names.begin(), it))] =
          options.fixed_parameter_values_m[lock];
    }

    TunaAssessmentControls controls = base_controls;
    bool start_ok = true;
    bool start_converged = true;
    int start_iters = 0;
    double start_gradient_norm = std::numeric_limits<double>::quiet_NaN();
    std::ostringstream start_msg;
    std::vector<TunaFitResult::PhaseDiagnostic> start_phase_diags;

    for (TunaAssessmentPhase phase : options.phase_sequence_m) {
      if (options.lbfgs_print_every_m > 0) {
        std::cout << "tuna_fit start=" << (start_idx + 1) << '/'
                  << options.multistart_m << " phase=" << phase_name(phase)
                  << " begin\n"
                  << std::flush;
      }
      controls = tuned_controls_for_phase(base_controls, phase);

      std::vector<bool> free_mask = build_phase_free_mask(fixed_names, phase);
      for (const std::string &locked_name : options.fixed_parameter_names_m) {
        const auto it =
            std::find(fixed_names.begin(), fixed_names.end(), locked_name);
        free_mask[static_cast<size_t>(std::distance(fixed_names.begin(), it))] =
            false;
      }
      const auto build_locked_from_mask = [&](const std::vector<bool> &mask,
                                              std::vector<size_t> &indices,
                                              std::vector<double> &values) {
        indices.clear();
        values.clear();
        indices.reserve(mask.size());
        values.reserve(mask.size());

        for (size_t j = 0; j < mask.size(); ++j) {
          if (!mask[j]) {
            indices.push_back(partition.fixed_indices_m[j]);
            values.push_back(theta[j]);
          }
        }
      };

      const auto absorb_fit = [&](const LaplaceExactLBFGSResult &fit,
                                  bool require_convergence) {
        ++best_fit.objective_evaluations_m;
        start_iters += fit.iterations_m;

        if (!fit.theta_hat_m.empty()) {
          theta = fit.theta_hat_m;
        }
        if (!fit.u_hat_m.empty()) {
          random = fit.u_hat_m;
        }

        start_gradient_norm = fit.gradient_norm_m;

        if (!std::isfinite(fit.laplace_objective_m)) {
          start_ok = false;
        }

        if (require_convergence && (!fit.converged_m || !fit.logdet_ok_m)) {
          start_converged = false;
        }
      };

      const auto absorb_joint_fit = [&](const JointOnlyExactLBFGSResult &fit,
                                        bool require_convergence) {
        ++best_fit.objective_evaluations_m;
        start_iters += fit.iterations_m;

        if (!fit.theta_hat_m.empty()) {
          theta = fit.theta_hat_m;
        }
        if (!fit.u_hat_m.empty()) {
          random = fit.u_hat_m;
        }

        start_gradient_norm = fit.gradient_norm_m;

        if (!std::isfinite(fit.joint_objective_m)) {
          start_ok = false;
        }

        if (require_convergence && !fit.converged_m) {
          start_converged = false;
        }
      };

      const auto record_diag = [&](const TunaAssessmentPhase diag_phase,
                                   const std::vector<bool> &diag_mask,
                                   const LaplaceExactLBFGSResult &fit,
                                   const std::string &tag = "") {
        TunaFitResult::PhaseDiagnostic diag;
        diag.phase_m = diag_phase;
        diag.start_index_m = start_idx;
        diag.free_fixed_count_m = static_cast<int>(
            std::count(diag_mask.begin(), diag_mask.end(), true));
        diag.locked_fixed_count_m =
            static_cast<int>(diag_mask.size()) - diag.free_fixed_count_m;
        diag.iterations_m = fit.iterations_m;
        diag.laplace_objective_m = fit.laplace_objective_m;
        diag.gradient_norm_m = fit.gradient_norm_m;
        std::vector<size_t> gradient_order(fit.gradient_fixed_m.size());
        std::iota(gradient_order.begin(), gradient_order.end(), size_t{0});
        std::sort(gradient_order.begin(), gradient_order.end(),
                  [&](size_t lhs, size_t rhs) {
                    return std::abs(fit.gradient_fixed_m[lhs]) >
                           std::abs(fit.gradient_fixed_m[rhs]);
                  });
        const size_t reported_gradient_count =
            std::min<size_t>(5, gradient_order.size());
        for (size_t rank = 0; rank < reported_gradient_count; ++rank) {
          const size_t j = gradient_order[rank];
          diag.largest_gradient_parameters_m.push_back(
              j < fixed_names.size() ? fixed_names[j]
                                     : "fixed_" + std::to_string(j));
          diag.largest_gradient_values_m.push_back(fit.gradient_fixed_m[j]);
        }
        diag.converged_m = fit.converged_m;
        diag.logdet_ok_m = fit.logdet_ok_m;
        if (tag.empty()) {
          diag.message_m = fit.message_m;
        } else {
          diag.message_m = tag + ": " + fit.message_m;
        }
        start_phase_diags.push_back(diag);

        start_msg << "[phase=" << static_cast<int>(diag_phase);
        if (!tag.empty()) {
          start_msg << " " << tag;
        }
        start_msg << " converged=" << (fit.converged_m ? "yes" : "no")
                  << " logdet=" << (fit.logdet_ok_m ? "ok" : "fail")
                  << " obj=" << fit.laplace_objective_m << "] ";
      };

      const auto record_joint_diag = [&](const TunaAssessmentPhase diag_phase,
                                         const std::vector<bool> &diag_mask,
                                         const JointOnlyExactLBFGSResult &fit,
                                         const std::string &tag = "") {
        TunaFitResult::PhaseDiagnostic diag;
        diag.phase_m = diag_phase;
        diag.start_index_m = start_idx;
        diag.free_fixed_count_m = static_cast<int>(
            std::count(diag_mask.begin(), diag_mask.end(), true));
        diag.locked_fixed_count_m =
            static_cast<int>(diag_mask.size()) - diag.free_fixed_count_m;
        diag.iterations_m = fit.iterations_m;
        diag.laplace_objective_m = fit.joint_objective_m;
        diag.gradient_norm_m = fit.gradient_norm_m;
        diag.converged_m = fit.converged_m;
        diag.logdet_ok_m = true;
        if (tag.empty()) {
          diag.message_m = "joint_only: " + fit.message_m;
        } else {
          diag.message_m = tag + ": " + fit.message_m;
        }
        start_phase_diags.push_back(diag);

        start_msg << "[phase=" << static_cast<int>(diag_phase);
        if (!tag.empty()) {
          start_msg << " " << tag;
        }
        start_msg << " converged=" << (fit.converged_m ? "yes" : "no")
                  << " logdet=na" << " obj=" << fit.joint_objective_m << "] ";
      };

      std::vector<size_t> locked_full_indices;
      std::vector<double> locked_full_values;
      build_locked_from_mask(free_mask, locked_full_indices,
                             locked_full_values);

      FixedTrustRegionSpec movement_entry_trust;
      if (phase == TunaAssessmentPhase::InitializeMovement) {
        movement_entry_trust = build_movement_entry_trust_region_spec(
            fixed_names, partition, theta, 0.35, 0.25);
      }

      PhaseLockedSpatialModel phase_model(
          data, controls, locked_full_indices, locked_full_values, {}, 1.0,
          movement_entry_trust.full_indices_m, movement_entry_trust.centers_m,
          movement_entry_trust.half_widths_m);
      const LaplaceExactLBFGSOptions lbfgs_options =
          make_exact_lbfgs_options(options);

      LaplaceExactLBFGSResult phase_fit;
      if (phase == TunaAssessmentPhase::Full) {
        phase_fit = optimize_full_phase_public_exact_laplace(
            phase_model, theta, random, partition, options, fixed_names);
      } else {
        phase_fit = optimize_laplace_fixed_effects_exact_lbfgs(
            phase_model, theta, random, partition, lbfgs_options);
      }

      const bool movement_init_failure =
          (phase == TunaAssessmentPhase::InitializeMovement) &&
          is_initial_exact_gradient_failure(phase_fit);

      // Initialization phases only provide warm starts.  Their stopping
      // status must not permanently determine whether the completed fit
      // converged; formal convergence is a property of the Full phase.
      absorb_fit(phase_fit,
                 phase == TunaAssessmentPhase::Full && !movement_init_failure);
      record_diag(phase, free_mask, phase_fit);
      if (options.lbfgs_print_every_m > 0) {
        std::cout << "tuna_fit start=" << (start_idx + 1) << '/'
                  << options.multistart_m << " phase=" << phase_name(phase)
                  << " end iterations=" << phase_fit.iterations_m
                  << " converged=" << (phase_fit.converged_m ? "true" : "false")
                  << '\n'
                  << std::flush;
      }

      if (movement_init_failure) {
        const std::vector<bool> bridge_mask = build_phase_free_mask(
            fixed_names, TunaAssessmentPhase::InitializeCatchability);
        std::vector<size_t> bridge_locked_indices;
        std::vector<double> bridge_locked_values;
        build_locked_from_mask(bridge_mask, bridge_locked_indices,
                               bridge_locked_values);

        TunaAssessmentControls bridge_controls = tuned_controls_for_phase(
            base_controls, TunaAssessmentPhase::InitializeCatchability);

        LaplaceExactLBFGSOptions bridge_options =
            make_exact_lbfgs_options(options);
        bridge_options.max_iterations_m =
            std::max(8, std::min(options.max_iterations_per_phase_m / 2, 20));
        bridge_options.initial_step_scale_m =
            std::min(bridge_options.initial_step_scale_m, 0.05);
        bridge_options.min_step_scale_m =
            std::min(bridge_options.min_step_scale_m, 1e-8);

        PhaseLockedSpatialModel bridge_model(
            data, bridge_controls, bridge_locked_indices, bridge_locked_values);

        LaplaceExactLBFGSResult bridge_fit =
            optimize_laplace_fixed_effects_exact_lbfgs(
                bridge_model, theta, random, partition, bridge_options);

        absorb_fit(bridge_fit, false);
        record_diag(phase, bridge_mask, bridge_fit, "bridge_lock");

        if (start_ok) {
          std::vector<double> movement_retry_random = random;
          {
            build_locked_from_mask(free_mask, locked_full_indices,
                                   locked_full_values);

            PhaseLockedSpatialModel warmup_model(
                data, controls, locked_full_indices, locked_full_values, {},
                1.0, movement_entry_trust.full_indices_m,
                movement_entry_trust.centers_m,
                movement_entry_trust.half_widths_m);

            const RandomEffectNewtonOptions warmup_options =
                make_movement_warmup_newton_options(options);

            RandomEffectNewtonResult warmup = optimize_random_effects_newton(
                warmup_model, theta, movement_retry_random, partition,
                warmup_options);

            ++best_fit.objective_evaluations_m;

            TunaFitResult::PhaseDiagnostic warmup_diag;
            warmup_diag.phase_m = phase;
            warmup_diag.start_index_m = start_idx;
            warmup_diag.free_fixed_count_m = static_cast<int>(
                std::count(free_mask.begin(), free_mask.end(), true));
            warmup_diag.locked_fixed_count_m =
                static_cast<int>(free_mask.size()) -
                warmup_diag.free_fixed_count_m;
            warmup_diag.iterations_m = warmup.iterations_m;
            warmup_diag.laplace_objective_m = warmup.objective_value_m;
            warmup_diag.gradient_norm_m = warmup.gradient_norm_m;
            warmup_diag.converged_m = warmup.converged_m;
            warmup_diag.logdet_ok_m = true;
            warmup_diag.message_m =
                "movement_warmup_newton: " + warmup.message_m;
            start_phase_diags.push_back(warmup_diag);

            start_msg << "[phase=" << static_cast<int>(phase)
                      << " movement_warmup_newton"
                      << " converged=" << (warmup.converged_m ? "yes" : "no")
                      << " logdet=na" << " obj=" << warmup.objective_value_m
                      << "] ";

            if (!warmup.u_hat_m.empty() &&
                warmup.u_hat_m.size() == movement_retry_random.size()) {
              movement_retry_random = warmup.u_hat_m;
              random = warmup.u_hat_m;
            }
          }

          const std::vector<bool> partial_move_mask =
              build_partial_movement_free_mask(fixed_names);
          build_locked_from_mask(partial_move_mask, locked_full_indices,
                                 locked_full_values);

          LaplaceExactLBFGSOptions partial_retry_options =
              make_exact_lbfgs_options(options);
          partial_retry_options.initial_step_scale_m =
              std::min(partial_retry_options.initial_step_scale_m, 0.03);
          partial_retry_options.min_step_scale_m =
              std::min(partial_retry_options.min_step_scale_m, 1e-8);
          partial_retry_options.max_iterations_m =
              std::max(10, std::min(options.max_iterations_per_phase_m, 30));

          PhaseLockedSpatialModel partial_retry_model(
              data, controls, locked_full_indices, locked_full_values, {}, 1.0,
              movement_entry_trust.full_indices_m,
              movement_entry_trust.centers_m,
              movement_entry_trust.half_widths_m);

          LaplaceExactLBFGSResult partial_retry_fit =
              optimize_laplace_fixed_effects_exact_lbfgs(
                  partial_retry_model, theta, movement_retry_random, partition,
                  partial_retry_options);

          absorb_fit(partial_retry_fit, false);
          record_diag(phase, partial_move_mask, partial_retry_fit,
                      "movement_partial_retry");

          const bool partial_ok =
              partial_retry_fit.logdet_ok_m &&
              !is_initial_exact_gradient_failure(partial_retry_fit);

          bool continuation_ok = partial_ok;
          bool used_joint_only_fallback = false;

          if (start_ok && !continuation_ok) {
            const std::vector<double> continuation_factors = {0.25, 0.5, 1.0};
            for (double factor : continuation_factors) {
              build_locked_from_mask(free_mask, locked_full_indices,
                                     locked_full_values);

              LaplaceExactLBFGSOptions ramp_options =
                  make_exact_lbfgs_options(options);
              ramp_options.initial_step_scale_m =
                  std::min(ramp_options.initial_step_scale_m, 0.03);
              ramp_options.min_step_scale_m =
                  std::min(ramp_options.min_step_scale_m, 1e-8);
              ramp_options.max_iterations_m =
                  std::max(8, std::min(options.max_iterations_per_phase_m, 20));

              PhaseLockedSpatialModel ramp_model(
                  data, controls, locked_full_indices, locked_full_values,
                  movement_availability_full_indices, factor,
                  movement_entry_trust.full_indices_m,
                  movement_entry_trust.centers_m,
                  movement_entry_trust.half_widths_m);

              LaplaceExactLBFGSResult ramp_fit =
                  optimize_laplace_fixed_effects_exact_lbfgs(
                      ramp_model, theta, movement_retry_random, partition,
                      ramp_options);

              absorb_fit(ramp_fit, false);

              const int factor_pct =
                  static_cast<int>(std::lround(factor * 100.0));
              record_diag(phase, free_mask, ramp_fit,
                          "movement_ramp_" + std::to_string(factor_pct));

              continuation_ok = ramp_fit.logdet_ok_m &&
                                !is_initial_exact_gradient_failure(ramp_fit);
              if (!start_ok || !continuation_ok) {
                if (!start_ok) {
                  break;
                }
              }
            }
          }

          if (start_ok && !partial_ok) {
            build_locked_from_mask(free_mask, locked_full_indices,
                                   locked_full_values);

            LaplaceExactLBFGSOptions reseed_options =
                make_exact_lbfgs_options(options);
            reseed_options.initial_step_scale_m =
                std::min(reseed_options.initial_step_scale_m, 0.03);
            reseed_options.min_step_scale_m =
                std::min(reseed_options.min_step_scale_m, 1e-8);
            reseed_options.max_iterations_m =
                std::max(10, std::min(options.max_iterations_per_phase_m, 35));
            reseed_options.warm_start_random_m = false;

            PhaseLockedSpatialModel reseed_model(
                data, controls, locked_full_indices, locked_full_values, {},
                1.0, movement_entry_trust.full_indices_m,
                movement_entry_trust.centers_m,
                movement_entry_trust.half_widths_m);

            LaplaceExactLBFGSResult reseed_fit =
                optimize_laplace_fixed_effects_exact_lbfgs(
                    reseed_model, theta, movement_retry_random, partition,
                    reseed_options);

            absorb_fit(reseed_fit, false);
            record_diag(phase, free_mask, reseed_fit, "movement_random_reseed");

            continuation_ok = continuation_ok ||
                              (reseed_fit.logdet_ok_m &&
                               !is_initial_exact_gradient_failure(reseed_fit));
          }

          if (start_ok && !continuation_ok) {
            build_locked_from_mask(free_mask, locked_full_indices,
                                   locked_full_values);

            JointOnlyExactLBFGSOptions joint_options =
                make_robust_joint_only_lbfgs_options(options);

            std::vector<double> random_zero(random.size(), 0.0);
            PhaseLockedSpatialModel joint_model(
                data, controls, locked_full_indices, locked_full_values, {},
                1.0, movement_entry_trust.full_indices_m,
                movement_entry_trust.centers_m,
                movement_entry_trust.half_widths_m);

            JointOnlyExactLBFGSResult joint_fit =
                optimize_joint_only_exact_lbfgs(joint_model, theta,
                                                movement_retry_random,
                                                partition, joint_options);

            absorb_joint_fit(joint_fit, false);
            record_joint_diag(phase, free_mask, joint_fit,
                              "movement_joint_stabilize");

            continuation_ok = !joint_fit.theta_hat_m.empty() &&
                              !joint_fit.u_hat_m.empty() &&
                              std::isfinite(joint_fit.joint_objective_m);
            used_joint_only_fallback = continuation_ok;
          }

          if (start_ok && continuation_ok && !used_joint_only_fallback) {
            build_locked_from_mask(free_mask, locked_full_indices,
                                   locked_full_values);

            LaplaceExactLBFGSOptions full_retry_options =
                make_exact_lbfgs_options(options);
            full_retry_options.initial_step_scale_m =
                std::min(full_retry_options.initial_step_scale_m, 0.05);
            full_retry_options.min_step_scale_m =
                std::min(full_retry_options.min_step_scale_m, 1e-8);

            PhaseLockedSpatialModel full_retry_model(
                data, controls, locked_full_indices, locked_full_values, {},
                1.0, movement_entry_trust.full_indices_m,
                movement_entry_trust.centers_m,
                movement_entry_trust.half_widths_m);

            LaplaceExactLBFGSResult full_retry_fit =
                optimize_laplace_fixed_effects_exact_lbfgs(
                    full_retry_model, theta, movement_retry_random, partition,
                    full_retry_options);

            absorb_fit(full_retry_fit, false);
            record_diag(phase, free_mask, full_retry_fit, "movement_retry");
          }
        }
      }

      if (!start_ok) {
        break;
      }
    }

    if (!start_ok) {
      continue;
    }

    const std::vector<double> full_hat =
        merge_parameters(theta, random, partition);

    TunaAssessmentControls final_controls = base_controls;
    final_controls.phase_m = TunaAssessmentPhase::Full;

    const TunaAssessmentRunSummary final_summary =
        evaluate_at_parameters(data, final_controls, full_hat,
                               "fitted_start_" + std::to_string(start_idx));

    const double final_obj = final_summary.nll_m;

    if (final_obj < best_obj) {
      best_obj = final_obj;
      best_fit.converged_m = start_converged;
      best_fit.best_start_index_m = start_idx;
      best_fit.total_iterations_m = start_iters;
      best_fit.gradient_norm_m = start_gradient_norm;
      best_fit.best_parameters_m = full_hat;
      best_fit.summary_m = final_summary;
      best_fit.message_m = start_msg.str();
      best_fit.phase_diagnostics_m = start_phase_diags;
    }
  }

  if (best_fit.best_start_index_m < 0) {
    best_fit.converged_m = false;
    best_fit.message_m =
        "fit_spatial_assessment: no successful start; using initial values";
    best_fit.summary_m =
        evaluate_initial_run(data, base_controls, "fitted_fallback_initial");
    best_fit.best_parameters_m = full_initial;
    best_fit.total_iterations_m = 0;
    best_fit.gradient_norm_m = std::numeric_limits<double>::quiet_NaN();
    best_fit.phase_diagnostics_m.clear();
  } else {
    best_fit.summary_m.label_m = "fitted";
  }

  return best_fit;
}

inline TunaSpatialAssessmentData
apply_sensitivity_scenario(const TunaSpatialAssessmentData &base,
                           const TunaSensitivityScenario &scenario) {
  TunaSpatialAssessmentData data = base;

  for (double &m : data.natural_mortality_at_age_m) {
    m *= scenario.natural_mortality_multiplier_m;
  }

  for (double &e : data.effort_m) {
    e *= scenario.effort_multiplier_m;
  }

  for (double &idx : data.observed_index_m) {
    idx *= scenario.index_multiplier_m;
  }

  for (double &a : data.availability_surface_m) {
    a *= scenario.availability_multiplier_m;
  }

  data.validate();
  return data;
}

inline std::vector<TunaAssessmentRunSummary>
run_sensitivity_initial(const TunaSpatialAssessmentData &base,
                        const TunaAssessmentControls &controls,
                        const std::vector<TunaSensitivityScenario> &scenarios) {
  std::vector<TunaAssessmentRunSummary> out;
  out.reserve(scenarios.size());

  for (const auto &scenario : scenarios) {
    const TunaSpatialAssessmentData modified =
        apply_sensitivity_scenario(base, scenario);
    out.push_back(evaluate_initial_run(modified, controls, scenario.label_m));
  }

  return out;
}

inline std::vector<TunaAssessmentRunSummary>
run_sensitivity_fit(const TunaSpatialAssessmentData &base,
                    const TunaAssessmentControls &controls,
                    const TunaFitOptions &fit_options,
                    const std::vector<TunaSensitivityScenario> &scenarios) {
  std::vector<TunaAssessmentRunSummary> out;
  out.reserve(scenarios.size());

  for (const auto &scenario : scenarios) {
    const TunaSpatialAssessmentData modified =
        apply_sensitivity_scenario(base, scenario);
    TunaFitResult fit = fit_spatial_assessment(modified, controls, fit_options);
    fit.summary_m.label_m = scenario.label_m;
    out.push_back(fit.summary_m);
  }

  return out;
}

inline TunaRetrospectiveResult
run_retrospective_initial(const TunaSpatialAssessmentData &base,
                          const TunaAssessmentControls &controls,
                          int max_peel) {
  if (max_peel < 0) {
    throw std::invalid_argument(
        "run_retrospective_initial: max_peel must be non-negative");
  }

  const int capped_peel = std::min(max_peel, std::max(0, base.n_years_m - 2));

  TunaRetrospectiveResult out;
  out.points_m.reserve(static_cast<size_t>(capped_peel + 1));

  for (int peel = 0; peel <= capped_peel; ++peel) {
    const int keep_years = base.n_years_m - peel;
    const TunaSpatialAssessmentData peeled = base.peel_years(keep_years);

    TunaRetrospectivePoint point;
    point.peel_m = peel;
    point.summary_m = evaluate_initial_run(
        peeled, controls, "retro_peel_" + std::to_string(peel));
    out.points_m.push_back(point);
  }

  if (!out.points_m.empty()) {
    const TunaAssessmentRunSummary &reference = out.points_m.front().summary_m;

    double rho_sum = 0.0;
    int rho_n = 0;
    for (size_t i = 1; i < out.points_m.size(); ++i) {
      const TunaAssessmentRunSummary &peeled = out.points_m[i].summary_m;
      const int ref_year_idx = peeled.n_years_m - 1;
      if (ref_year_idx < 0 ||
          ref_year_idx >= static_cast<int>(reference.ssb_by_year_m.size())) {
        continue;
      }

      const double ref_ssb =
          reference.ssb_by_year_m[static_cast<size_t>(ref_year_idx)];
      const double peel_ssb = peeled.ssb_terminal_m;
      if (std::isfinite(ref_ssb) && std::fabs(ref_ssb) > 1e-12 &&
          std::isfinite(peel_ssb)) {
        rho_sum += (peel_ssb - ref_ssb) / ref_ssb;
        ++rho_n;
      }
    }

    if (rho_n > 0) {
      out.mohns_rho_ssb_m = rho_sum / static_cast<double>(rho_n);
    }
  }

  return out;
}

inline TunaRetrospectiveResult
run_retrospective_fit(const TunaSpatialAssessmentData &base,
                      const TunaAssessmentControls &controls,
                      const TunaFitOptions &fit_options, int max_peel) {
  if (max_peel < 0) {
    throw std::invalid_argument(
        "run_retrospective_fit: max_peel must be non-negative");
  }

  const int capped_peel = std::min(max_peel, std::max(0, base.n_years_m - 2));

  TunaRetrospectiveResult out;
  out.points_m.reserve(static_cast<size_t>(capped_peel + 1));

  for (int peel = 0; peel <= capped_peel; ++peel) {
    const int keep_years = base.n_years_m - peel;
    const TunaSpatialAssessmentData peeled = base.peel_years(keep_years);

    TunaRetrospectivePoint point;
    point.peel_m = peel;
    TunaFitResult fit = fit_spatial_assessment(peeled, controls, fit_options);
    fit.summary_m.label_m = "retro_fit_peel_" + std::to_string(peel);
    point.summary_m = fit.summary_m;
    out.points_m.push_back(point);
  }

  if (!out.points_m.empty()) {
    const TunaAssessmentRunSummary &reference = out.points_m.front().summary_m;
    double rho_sum = 0.0;
    int rho_n = 0;

    for (size_t i = 1; i < out.points_m.size(); ++i) {
      const TunaAssessmentRunSummary &peeled = out.points_m[i].summary_m;
      const int ref_year_idx = peeled.n_years_m - 1;
      if (ref_year_idx < 0 ||
          ref_year_idx >= static_cast<int>(reference.ssb_by_year_m.size())) {
        continue;
      }

      const double ref_ssb =
          reference.ssb_by_year_m[static_cast<size_t>(ref_year_idx)];
      const double peel_ssb = peeled.ssb_terminal_m;
      if (std::isfinite(ref_ssb) && std::fabs(ref_ssb) > 1e-12 &&
          std::isfinite(peel_ssb)) {
        rho_sum += (peel_ssb - ref_ssb) / ref_ssb;
        ++rho_n;
      }
    }

    if (rho_n > 0) {
      out.mohns_rho_ssb_m = rho_sum / static_cast<double>(rho_n);
    }
  }

  return out;
}

inline TunaSpatialAssessmentData
simulate_observed_data(const TunaSpatialAssessmentData &truth,
                       std::uint64_t seed, double index_cv = 0.2,
                       double biomass_cv = 0.2) {
  TunaSpatialAssessmentData sim = truth;
  std::mt19937_64 rng(seed);

  std::normal_distribution<double> z(0.0, 1.0);
  for (double &idx : sim.observed_index_m) {
    const double sigma = std::max(1e-6, index_cv);
    const double eps = z(rng);
    idx = std::max(1e-6, idx * std::exp(-0.5 * sigma * sigma + sigma * eps));
  }

  for (double &ret : sim.observed_retained_biomass_m) {
    const double sigma = std::max(1e-6, biomass_cv);
    const double eps = z(rng);
    ret = std::max(1e-6, ret * std::exp(-0.5 * sigma * sigma + sigma * eps));
  }

  for (double &dis : sim.observed_discard_biomass_m) {
    const double sigma = std::max(1e-6, biomass_cv);
    const double eps = z(rng);
    dis = std::max(1e-6, dis * std::exp(-0.5 * sigma * sigma + sigma * eps));
  }

  std::poisson_distribution<int> pois;
  for (int &c : sim.observed_catch_numbers_m) {
    const double lambda = std::max(1.0, static_cast<double>(c));
    pois = std::poisson_distribution<int>(lambda);
    c = std::max(0, pois(rng));
  }

  sim.validate();
  return sim;
}

inline TunaSimulationResult
run_simulation_estimation_loop(const TunaSpatialAssessmentData &truth,
                               const TunaAssessmentControls &controls,
                               const TunaFitOptions &fit_options,
                               int n_simulations, std::uint64_t seed) {
  if (n_simulations <= 0) {
    throw std::invalid_argument(
        "run_simulation_estimation_loop: n_simulations must be positive");
  }

  TunaSimulationResult out;
  out.cases_m.reserve(static_cast<size_t>(n_simulations));
  std::vector<double> finite_biases;
  finite_biases.reserve(static_cast<size_t>(n_simulations));
  int low_depletion_count = 0;

  double bias_sum = 0.0;
  int bias_n = 0;

  for (int sim_id = 0; sim_id < n_simulations; ++sim_id) {
    TunaSimulationCase sim_case;
    sim_case.simulation_id_m = sim_id + 1;

    const TunaSpatialAssessmentData synthetic = simulate_observed_data(
        truth, seed + static_cast<std::uint64_t>(sim_id));

    sim_case.truth_m = evaluate_initial_run(truth, controls, "truth_initial");

    TunaFitResult fit =
        fit_spatial_assessment(synthetic, controls, fit_options);
    fit.summary_m.label_m = "sim_" + std::to_string(sim_case.simulation_id_m);
    sim_case.estimated_m = fit.summary_m;

    if (std::isfinite(sim_case.truth_m.depletion_terminal_m) &&
        std::isfinite(sim_case.estimated_m.depletion_terminal_m)) {
      sim_case.depletion_bias_m = sim_case.estimated_m.depletion_terminal_m -
                                  sim_case.truth_m.depletion_terminal_m;
      bias_sum += sim_case.depletion_bias_m;
      ++bias_n;
      finite_biases.push_back(sim_case.depletion_bias_m);

      if (sim_case.estimated_m.depletion_terminal_m <
          out.low_depletion_threshold_m) {
        ++low_depletion_count;
      }
    }

    out.cases_m.push_back(sim_case);
  }

  if (bias_n > 0) {
    out.mean_depletion_bias_m = bias_sum / static_cast<double>(bias_n);
    out.n_finite_bias_m = bias_n;

    std::sort(finite_biases.begin(), finite_biases.end());
    out.median_depletion_bias_m = empirical_quantile(finite_biases, 0.50);
    out.p10_depletion_bias_m = empirical_quantile(finite_biases, 0.10);
    out.p90_depletion_bias_m = empirical_quantile(finite_biases, 0.90);
    out.low_depletion_count_m = low_depletion_count;
    out.low_depletion_rate_m =
        static_cast<double>(low_depletion_count) / static_cast<double>(bias_n);
  }

  return out;
}

inline std::string biomass_trajectory_csv(const TunaAssessmentRunSummary &run) {
  std::ostringstream ss;
  ss << "year,spawning_biomass,depletion\n";
  for (size_t y = 0; y < run.ssb_by_year_m.size(); ++y) {
    const double ssb = run.ssb_by_year_m[y];
    const double depletion = std::isfinite(run.ssb0_m) && run.ssb0_m > 0.0
                                 ? ssb / run.ssb0_m
                                 : std::numeric_limits<double>::quiet_NaN();
    ss << (y + 1) << "," << ssb << "," << depletion << "\n";
  }
  return ss.str();
}

inline std::string
spatial_animation_csv(const TunaSpatialAssessmentData &data,
                      const TunaAssessmentControls &base_controls,
                      const std::vector<double> &parameters) {
  TunaAssessmentControls controls = base_controls;
  controls.phase_m = TunaAssessmentPhase::Full;
  controls.report_spatial_animation_m = true;
  AdvancedSpatialTunaAssessmentModel model(data, controls);
  if (parameters.size() != model.parameter_set().size())
    throw std::invalid_argument(
        "spatial_animation_csv: parameter length mismatch");

  ModelReportContext ctx;
  model.evaluate(parameters, ctx);
  std::ostringstream ss;
  ss << "record_type,year,season,region,fleet,biomass,retained_catch,"
        "discard_catch,from_region,to_region,movement_biomass,"
        "movement_probability\n";
  for (int y = 0; y < data.n_years_m; ++y) {
    for (int s = 0; s < data.n_seasons_m; ++s) {
      const std::string frame =
          "y" + std::to_string(y + 1) + "_s" + std::to_string(s + 1);
      for (int r = 0; r < data.n_regions_m; ++r) {
        const double biomass = report_value_or_nan(
            ctx, "anim_biomass_" + frame + "_r" + std::to_string(r + 1));
        for (int f = 0; f < data.n_fleets_m; ++f) {
          const std::string fleet_region = frame + "_f" +
                                           std::to_string(f + 1) + "_r" +
                                           std::to_string(r + 1);
          ss << "region," << (y + 1) << "," << (s + 1) << "," << (r + 1) << ","
             << (f + 1) << "," << biomass << ","
             << report_value_or_nan(ctx, "anim_retained_" + fleet_region) << ","
             << report_value_or_nan(ctx, "anim_discard_" + fleet_region)
             << ",,,,\n";
        }
        for (int to = 0; to < data.n_regions_m; ++to) {
          const std::string route = frame + "_from_" + std::to_string(r + 1) +
                                    "_to_" + std::to_string(to + 1);
          const std::string probability =
              "move_prob_s" + std::to_string(s + 1) + "_from_" +
              std::to_string(r + 1) + "_to_" + std::to_string(to + 1);
          ss << "movement," << (y + 1) << "," << (s + 1) << ",,,,,," << (r + 1)
             << "," << (to + 1) << ","
             << report_value_or_nan(ctx, "anim_movement_" + route) << ","
             << report_value_or_nan(ctx, probability) << "\n";
        }
      }
    }
  }
  return ss.str();
}

inline std::string
management_summary_csv(const std::vector<TunaAssessmentRunSummary> &runs) {
  std::ostringstream ss;
  ss << "label,years,nll,r0,steepness,ssb0,ssb_terminal,depletion_terminal\n";
  for (const auto &run : runs) {
    ss << run.label_m << "," << run.n_years_m << "," << run.nll_m << ","
       << run.r0_m << "," << run.steepness_m << "," << run.ssb0_m << ","
       << run.ssb_terminal_m << "," << run.depletion_terminal_m << "\n";
  }
  return ss.str();
}

inline std::string
observation_diagnostics_csv(const TunaDiagnosticBundle &diagnostics) {
  std::ostringstream ss;
  ss << "component,year,season,fleet,region,age,observed,predicted,"
        "standard_deviation,standardized_residual\n";
  for (const auto &row : diagnostics.observations_m) {
    ss << row.component_m << "," << row.year_m << "," << row.season_m << ","
       << row.fleet_m << "," << row.region_m << "," << row.age_m << ","
       << row.observed_m << "," << row.predicted_m << ","
       << row.standard_deviation_m << "," << row.standardized_residual_m
       << "\n";
  }
  return ss.str();
}

inline std::string
residual_summary_csv(const TunaDiagnosticBundle &diagnostics) {
  std::ostringstream ss;
  ss << "component,fleet,n,mean_residual,sdnr,rmse\n";
  for (const auto &row : diagnostics.residual_summaries_m) {
    ss << row.component_m << "," << row.fleet_m << "," << row.n_m << ","
       << row.mean_residual_m << "," << row.sdnr_m << "," << row.rmse_m << "\n";
  }
  return ss.str();
}

inline std::string
stratified_residual_summary_csv(const TunaDiagnosticBundle &diagnostics) {
  std::ostringstream ss;
  ss << "component,fleet,stratum,level,n,mean_residual,sdnr,rmse\n";
  for (const auto &row : diagnostics.stratified_residuals_m) {
    ss << row.component_m << "," << row.fleet_m << "," << row.stratum_m << ","
       << row.level_m << "," << row.n_m << "," << row.mean_residual_m << ","
       << row.sdnr_m << "," << row.rmse_m << "\n";
  }
  return ss.str();
}

inline std::string
likelihood_decomposition_csv(const TunaDiagnosticBundle &diagnostics) {
  std::ostringstream ss;
  ss << "component,nll,absolute_share,observation_units,nll_per_unit\n";
  for (const auto &row : diagnostics.likelihood_components_m)
    ss << row.component_m << "," << row.nll_m << "," << row.share_absolute_m
       << "," << row.observation_units_m << "," << row.nll_per_unit_m << "\n";
  return ss.str();
}

inline std::string
parameter_diagnostics_csv(const TunaDiagnosticBundle &diagnostics) {
  std::ostringstream ss;
  ss << "name,value,initial,displacement,random_effect,extreme\n";
  for (const auto &row : diagnostics.parameters_m)
    ss << row.name_m << "," << row.value_m << "," << row.initial_m << ","
       << row.displacement_m << "," << (row.random_effect_m ? 1 : 0) << ","
       << (row.extreme_m ? 1 : 0) << "\n";
  return ss.str();
}

inline std::string
catchability_availability_csv(const TunaDiagnosticBundle &diagnostics) {
  std::ostringstream ss;
  ss << "fleet,season,region,log_fishing_q,log_index_q,log_availability,"
        "log_fishing_q_availability\n";
  for (const auto &row : diagnostics.catchability_availability_m)
    ss << row.fleet_m << "," << row.season_m << "," << row.region_m << ","
       << row.log_fishing_q_m << "," << row.log_index_q_m << ","
       << row.log_availability_m << "," << row.log_fishing_q_availability_m
       << "\n";
  return ss.str();
}

inline const TunaResidualSummary *
find_residual_summary(const TunaDiagnosticBundle &diagnostics,
                      const std::string &component, int fleet) {
  for (const auto &row : diagnostics.residual_summaries_m)
    if (row.component_m == component && row.fleet_m == fleet)
      return &row;
  return nullptr;
}

inline std::string composition_weight_grid_csv(
    const std::vector<double> &weights, const std::vector<TunaFitResult> &fits,
    const std::vector<TunaDiagnosticBundle> &diagnostics, int n_fleets) {
  if (weights.size() != fits.size() || weights.size() != diagnostics.size())
    throw std::invalid_argument(
        "composition_weight_grid_csv: grid vector length mismatch");

  std::ostringstream ss;
  ss << "weight,converged,nll,ssb0,ssb_terminal,depletion_terminal";
  for (int fleet = 1; fleet <= n_fleets; ++fleet)
    ss << ",index_sdnr_f" << fleet << ",index_rmse_f" << fleet
       << ",composition_sdnr_f" << fleet << ",composition_rmse_f" << fleet;
  ss << "\n";

  for (size_t i = 0; i < weights.size(); ++i) {
    const auto &fit = fits[i];
    ss << weights[i] << "," << (fit.converged_m ? 1 : 0) << ","
       << fit.summary_m.nll_m << "," << fit.summary_m.ssb0_m << ","
       << fit.summary_m.ssb_terminal_m << ","
       << fit.summary_m.depletion_terminal_m;
    for (int fleet = 1; fleet <= n_fleets; ++fleet) {
      const TunaResidualSummary *index =
          find_residual_summary(diagnostics[i], "index", fleet);
      const TunaResidualSummary *composition =
          find_residual_summary(diagnostics[i], "composition", fleet);
      ss << ","
         << (index ? index->sdnr_m : std::numeric_limits<double>::quiet_NaN())
         << ","
         << (index ? index->rmse_m : std::numeric_limits<double>::quiet_NaN())
         << ","
         << (composition ? composition->sdnr_m
                         : std::numeric_limits<double>::quiet_NaN())
         << ","
         << (composition ? composition->rmse_m
                         : std::numeric_limits<double>::quiet_NaN());
    }
    ss << "\n";
  }
  return ss.str();
}

inline std::string
retrospective_summary_csv(const TunaRetrospectiveResult &retro) {
  std::ostringstream ss;
  ss << "peel,years,nll,depletion_terminal\n";
  for (const auto &point : retro.points_m) {
    ss << point.peel_m << "," << point.summary_m.n_years_m << ","
       << point.summary_m.nll_m << "," << point.summary_m.depletion_terminal_m
       << "\n";
  }
  ss << "mohns_rho_ssb,,," << retro.mohns_rho_ssb_m << "\n";
  return ss.str();
}

inline std::string simulation_summary_csv(const TunaSimulationResult &sim) {
  std::ostringstream ss;
  ss << "simulation_id,truth_depletion,estimated_depletion,depletion_bias\n";
  for (const auto &c : sim.cases_m) {
    ss << c.simulation_id_m << "," << c.truth_m.depletion_terminal_m << ","
       << c.estimated_m.depletion_terminal_m << "," << c.depletion_bias_m
       << "\n";
  }
  ss << "mean_bias,,," << sim.mean_depletion_bias_m << "\n";
  ss << "median_bias,,," << sim.median_depletion_bias_m << "\n";
  ss << "p10_bias,,," << sim.p10_depletion_bias_m << "\n";
  ss << "p90_bias,,," << sim.p90_depletion_bias_m << "\n";
  ss << "low_depletion_rate,,," << sim.low_depletion_rate_m << "\n";
  ss << "low_depletion_threshold,,," << sim.low_depletion_threshold_m << "\n";
  ss << "low_depletion_count,,," << sim.low_depletion_count_m << "\n";
  ss << "n_finite_bias,,," << sim.n_finite_bias_m << "\n";
  return ss.str();
}

inline std::string json_escape(const std::string &s) {
  std::ostringstream out;
  for (char c : s) {
    switch (c) {
    case '\\':
      out << "\\\\";
      break;
    case '"':
      out << "\\\"";
      break;
    case '\n':
      out << "\\n";
      break;
    case '\r':
      out << "\\r";
      break;
    case '\t':
      out << "\\t";
      break;
    default:
      out << c;
      break;
    }
  }
  return out.str();
}

inline std::string json_number(double value) {
  if (!std::isfinite(value)) {
    return "null";
  }
  std::ostringstream out;
  out << std::setprecision(17) << value;
  return out.str();
}

inline std::string
fit_result_json(const TunaFitResult &fit, const TunaRetrospectiveResult &retro,
                const std::vector<TunaAssessmentRunSummary> &sensitivity) {
  std::ostringstream ss;
  ss << "{\n";
  ss << "  \"fit\": {\n";
  ss << "    \"converged\": " << (fit.converged_m ? "true" : "false") << ",\n";
  ss << "    \"best_start_index\": " << fit.best_start_index_m << ",\n";
  ss << "    \"objective_evaluations\": " << fit.objective_evaluations_m
     << ",\n";
  ss << "    \"total_iterations\": " << fit.total_iterations_m << ",\n";
  ss << "    \"gradient_norm\": " << json_number(fit.gradient_norm_m) << ",\n";
  ss << "    \"message\": \"" << json_escape(fit.message_m) << "\",\n";
  ss << "    \"nll\": " << json_number(fit.summary_m.nll_m) << ",\n";
  ss << "    \"depletion_terminal\": "
     << json_number(fit.summary_m.depletion_terminal_m) << "\n";
  ss << "  },\n";

  ss << "  \"phase_diagnostics\": [\n";
  for (size_t i = 0; i < fit.phase_diagnostics_m.size(); ++i) {
    const auto &d = fit.phase_diagnostics_m[i];
    ss << "    {\"phase\": \"" << phase_name(d.phase_m)
       << "\", \"start_index\": " << d.start_index_m
       << ", \"free_fixed_count\": " << d.free_fixed_count_m
       << ", \"locked_fixed_count\": " << d.locked_fixed_count_m
       << ", \"iterations\": " << d.iterations_m
       << ", \"laplace_objective\": " << json_number(d.laplace_objective_m)
       << ", \"gradient_norm\": " << json_number(d.gradient_norm_m)
       << ", \"largest_gradients\": [";
    for (size_t j = 0; j < d.largest_gradient_parameters_m.size(); ++j) {
      ss << "{\"parameter\": \""
         << json_escape(d.largest_gradient_parameters_m[j])
         << "\", \"value\": " << json_number(d.largest_gradient_values_m[j])
         << "}";
      if (j + 1 < d.largest_gradient_parameters_m.size()) {
        ss << ", ";
      }
    }
    ss << "]" << ", \"converged\": " << (d.converged_m ? "true" : "false")
       << ", \"logdet_ok\": " << (d.logdet_ok_m ? "true" : "false")
       << ", \"message\": \"" << json_escape(d.message_m) << "\"}";
    if (i + 1 < fit.phase_diagnostics_m.size()) {
      ss << ",";
    }
    ss << "\n";
  }
  ss << "  ],\n";

  ss << "  \"retrospective\": {\n";
  ss << "    \"mohns_rho_ssb\": " << json_number(retro.mohns_rho_ssb_m)
     << ",\n";
  ss << "    \"points\": [\n";
  for (size_t i = 0; i < retro.points_m.size(); ++i) {
    const auto &p = retro.points_m[i];
    ss << "      {\"peel\": " << p.peel_m
       << ", \"years\": " << p.summary_m.n_years_m
       << ", \"depletion_terminal\": "
       << json_number(p.summary_m.depletion_terminal_m) << "}";
    if (i + 1 < retro.points_m.size()) {
      ss << ",";
    }
    ss << "\n";
  }
  ss << "    ]\n";
  ss << "  },\n";

  ss << "  \"sensitivity\": [\n";
  for (size_t i = 0; i < sensitivity.size(); ++i) {
    const auto &s = sensitivity[i];
    ss << "    {\"label\": \"" << json_escape(s.label_m)
       << "\", \"nll\": " << json_number(s.nll_m)
       << ", \"depletion_terminal\": " << json_number(s.depletion_terminal_m)
       << "}";
    if (i + 1 < sensitivity.size()) {
      ss << ",";
    }
    ss << "\n";
  }
  ss << "  ]\n";
  ss << "}\n";
  return ss.str();
}

inline void write_text_file(const std::string &path, const std::string &text) {
  std::ofstream out(path);
  if (!out) {
    throw std::runtime_error("write_text_file: unable to open " + path);
  }
  out << text;
}

} // namespace quadra

#endif // QUADRA_TUNA_ASSESSMENT_ACCEPTANCE_HPP
