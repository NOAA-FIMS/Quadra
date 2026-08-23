#include <algorithm>
#include <atomic>
#include <cctype>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <random>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include "../include/tuna/dependency_free_transport_flow.hpp"
#include <Eigen/Dense>

#include "../include/tuna/tuna_assessment_acceptance.hpp"
#include "../include/tuna/tuna_reference_points.hpp"
#include "../include/tuna/tuna_spatial_assessment_model.hpp"
#include "include/quadra/sampling.hpp"

DECLARE_ADGRAPH();

namespace {
class AssessmentConfig {
public:
  void load(const std::string &path) {
    std::ifstream input(path);
    if (!input)
      throw std::runtime_error("could not open assessment config: " + path);
    source_m = path;
    std::string line;
    int line_number = 0;
    while (std::getline(input, line)) {
      ++line_number;
      const size_t comment = line.find('#');
      if (comment != std::string::npos)
        line.erase(comment);
      line = trim(line);
      if (line.empty() || (line.front() == '[' && line.back() == ']'))
        continue;
      const size_t equals = line.find('=');
      const std::string key = trim(line.substr(0, equals));
      if (key.empty())
        continue;
      if (equals == std::string::npos)
        throw std::runtime_error(path + ":" + std::to_string(line_number) +
                                 ": expected key = value");
      const std::string value = trim(line.substr(equals + 1));
      if (value.empty())
        throw std::runtime_error(path + ":" + std::to_string(line_number) +
                                 ": configuration value is empty");
      if (!values_m.emplace(key, value).second)
        throw std::runtime_error(path + ":" + std::to_string(line_number) +
                                 ": duplicate configuration key " + key);
    }
  }

  std::string get(const std::string &key, const std::string &fallback) const {
    const auto found = values_m.find(key);
    return found == values_m.end() ? fallback : found->second;
  }

  bool has(const std::string &key) const {
    return values_m.find(key) != values_m.end();
  }

  const std::string &source() const { return source_m; }

private:
  static std::string trim(std::string value) {
    const auto nonspace = [](unsigned char c) { return !std::isspace(c); };
    value.erase(value.begin(),
                std::find_if(value.begin(), value.end(), nonspace));
    value.erase(std::find_if(value.rbegin(), value.rend(), nonspace).base(),
                value.end());
    return value;
  }

  std::unordered_map<std::string, std::string> values_m;
  std::string source_m;
};

AssessmentConfig g_config;

std::string configured_string(const char *environment, const std::string &key,
                              const std::string &fallback) {
  if (const char *value = std::getenv(environment))
    return value;
  return g_config.get(key, fallback);
}

int configured_int(const char *environment, const std::string &key,
                   int fallback) {
  return std::stoi(
      configured_string(environment, key, std::to_string(fallback)));
}

double configured_double(const char *environment, const std::string &key,
                         double fallback) {
  std::ostringstream text;
  text << std::setprecision(17) << fallback;
  return std::stod(configured_string(environment, key, text.str()));
}

bool configured_bool(const char *environment, const std::string &key,
                     bool fallback) {
  std::string value =
      configured_string(environment, key, fallback ? "true" : "false");
  std::transform(value.begin(), value.end(), value.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  if (value == "1" || value == "true" || value == "yes" || value == "on")
    return true;
  if (value == "0" || value == "false" || value == "no" || value == "off")
    return false;
  throw std::invalid_argument("invalid boolean for " + key + ": " + value);
}

void write_effective_configuration(const std::string &out_dir) {
  const std::vector<std::pair<std::string, std::string>> values = {
      {"config.source", g_config.source()},
      {"output.data_dir",
       configured_string("QUADRA_TUNA_OUTPUT_DIR", "output.data_dir", out_dir)},
      {"data.model_consistent",
       configured_bool("QUADRA_TUNA_MODEL_CONSISTENT_DATA",
                       "data.model_consistent", true)
           ? "true"
           : "false"},
      {"model.anchor_fleet",
       std::to_string(configured_int("QUADRA_TUNA_ANCHOR_FLEET",
                                     "model.anchor_fleet", 1))},
      {"fit.multistart", std::to_string(configured_int("QUADRA_TUNA_MULTISTART",
                                                       "fit.multistart", 4))},
      {"fit.max_phase_iterations",
       std::to_string(configured_int("QUADRA_TUNA_MAX_PHASE_ITERATIONS",
                                     "fit.max_phase_iterations", 35))},
      {"fit.hdot_workers",
       std::to_string(
           configured_int("QUADRA_TUNA_HDOT_WORKERS", "fit.hdot_workers", 0))},
      {"fit.load_checkpoint", configured_bool("QUADRA_TUNA_LOAD_FIT_CHECKPOINT",
                                              "fit.load_checkpoint", false)
                                  ? "true"
                                  : "false"},
      {"sampling.enabled",
       configured_bool("QUADRA_TUNA_RUN_NUTS", "sampling.enabled", false)
           ? "true"
           : "false"},
      {"sampling.method",
       configured_string("QUADRA_TUNA_SAMPLER", "sampling.method", "pcn")},
      {"sampling.chains",
       std::to_string(
           configured_int("QUADRA_TUNA_NUTS_CHAINS", "sampling.chains", 4))},
      {"sampling.warmup",
       std::to_string(
           configured_int("QUADRA_TUNA_NUTS_WARMUP", "sampling.warmup", 500))},
      {"sampling.samples",
       std::to_string(configured_int("QUADRA_TUNA_NUTS_SAMPLES",
                                     "sampling.samples", 500))},
      {"sampling.max_tree_depth",
       std::to_string(configured_int("QUADRA_TUNA_NUTS_MAX_TREE_DEPTH",
                                     "sampling.max_tree_depth", 10))},
      {"sampling.target_acceptance",
       std::to_string(configured_double("QUADRA_TUNA_NUTS_TARGET_ACCEPTANCE",
                                        "sampling.target_acceptance", 0.8))},
      {"sampling.initial_step_size",
       std::to_string(configured_double("QUADRA_TUNA_NUTS_INITIAL_STEP_SIZE",
                                        "sampling.initial_step_size", 0.0))},
      {"sampling.adapt_mass", configured_bool("QUADRA_TUNA_NUTS_ADAPT_MASS",
                                              "sampling.adapt_mass", true)
                                  ? "true"
                                  : "false"},
      {"sampling.dense_metric", configured_bool("QUADRA_TUNA_NUTS_DENSE_METRIC",
                                                "sampling.dense_metric", true)
                                    ? "true"
                                    : "false"},
      {"sampling.gradient_workers",
       std::to_string(configured_int("QUADRA_TUNA_GRADIENT_WORKERS",
                                     "sampling.gradient_workers", 4))},
      {"sampling.newton_initial_step",
       std::to_string(configured_double("QUADRA_TUNA_NUTS_NEWTON_INITIAL_STEP",
                                        "sampling.newton_initial_step", 1.0))},
      {"sampling.parallel_chains",
       configured_bool("QUADRA_TUNA_NUTS_PARALLEL_CHAINS",
                       "sampling.parallel_chains", true)
           ? "true"
           : "false"},
      {"sampling.tape_rebuild_interval",
       std::to_string(configured_int("QUADRA_TUNA_NUTS_TAPE_REBUILD_INTERVAL",
                                     "sampling.tape_rebuild_interval", 1))},
      {"sampling.pcn_beta",
       std::to_string(configured_double("QUADRA_TUNA_PCN_BETA",
                                        "sampling.pcn_beta", 0.5))},
      {"sampling.pcn_target_acceptance",
       std::to_string(configured_double("QUADRA_TUNA_PCN_TARGET_ACCEPTANCE",
                                        "sampling.pcn_target_acceptance",
                                        0.3))},
      {"sampling.pcn_adapt_beta",
       configured_bool("QUADRA_TUNA_PCN_ADAPT_BETA", "sampling.pcn_adapt_beta",
                       false)
           ? "true"
           : "false"},
      {"sampling.proposal_draws",
       configured_string(
           "QUADRA_TUNA_PROPOSAL_DRAWS", "sampling.proposal_draws",
           "build/assessment_outputs/benchmark_data/posterior_draws.csv")},
      {"sampling.proposal_df",
       std::to_string(configured_double("QUADRA_TUNA_PROPOSAL_DF",
                                        "sampling.proposal_df", 8.0))},
      {"sampling.proposal_inflation",
       std::to_string(configured_double("QUADRA_TUNA_PROPOSAL_INFLATION",
                                        "sampling.proposal_inflation", 1.2))},
      {"sampling.proposal_global_weight",
       std::to_string(configured_double("QUADRA_TUNA_PROPOSAL_GLOBAL_WEIGHT",
                                        "sampling.proposal_global_weight",
                                        0.9))},
      {"sampling.proposal_local_scale",
       std::to_string(configured_double("QUADRA_TUNA_PROPOSAL_LOCAL_SCALE",
                                        "sampling.proposal_local_scale", 0.2))},
      {"sampling.proposal_rw_scale",
       std::to_string(configured_double("QUADRA_TUNA_PROPOSAL_RW_SCALE",
                                        "sampling.proposal_rw_scale", 0.5))},
      {"sampling.transport_components",
       std::to_string(configured_int("QUADRA_TUNA_TRANSPORT_COMPONENTS",
                                     "sampling.transport_components", 6))},
      {"sampling.transport_inflation",
       std::to_string(configured_double("QUADRA_TUNA_TRANSPORT_INFLATION",
                                        "sampling.transport_inflation", 1.1))},
      {"sampling.transport_shrinkage",
       std::to_string(configured_double("QUADRA_TUNA_TRANSPORT_SHRINKAGE",
                                        "sampling.transport_shrinkage", 0.15))},
      {"sampling.transport_bandwidth",
       std::to_string(configured_double("QUADRA_TUNA_TRANSPORT_BANDWIDTH",
                                        "sampling.transport_bandwidth", 0.15))},
      {"sampling.transport_ridge",
       std::to_string(configured_double("QUADRA_TUNA_TRANSPORT_RIDGE",
                                        "sampling.transport_ridge", 0.001))},
      {"sampling.transport_model",
       configured_string("QUADRA_TUNA_TRANSPORT_MODEL",
                         "sampling.transport_model",
                         "build/transport_flow/native_realnvp_20260823.qflow")},
      {"sampling.transport_models",
       configured_string("QUADRA_TUNA_TRANSPORT_MODELS",
                         "sampling.transport_models",
                         "build/transport_flow/native_realnvp_20260823.qflow")},
      {"sampling.transport_candidates",
       std::to_string(configured_int("QUADRA_TUNA_TRANSPORT_CANDIDATES",
                                     "sampling.transport_candidates", 4))}};
  std::ostringstream csv;
  csv << "key,value\n";
  for (const auto &entry : values)
    csv << entry.first << "," << entry.second << "\n";
  quadra::write_text_file(out_dir + "/effective_configuration.csv", csv.str());
}

template <class Type>
void fingerprint_value(std::uint64_t &hash, const Type &value) {
  const unsigned char *bytes = reinterpret_cast<const unsigned char *>(&value);
  for (size_t i = 0; i < sizeof(Type); ++i) {
    hash ^= bytes[i];
    hash *= 1099511628211ULL;
  }
}

template <class Type>
void fingerprint_vector(std::uint64_t &hash, const std::vector<Type> &values) {
  fingerprint_value(hash, values.size());
  for (const Type &value : values)
    fingerprint_value(hash, value);
}

std::string
assessment_fingerprint(const quadra::TunaSpatialAssessmentData &data,
                       const quadra::TunaAssessmentControls &controls) {
  std::uint64_t hash = 1469598103934665603ULL;
  fingerprint_value(hash, data.n_years_m);
  fingerprint_value(hash, data.n_ages_m);
  fingerprint_value(hash, data.n_fleets_m);
  fingerprint_value(hash, data.n_regions_m);
  fingerprint_value(hash, data.n_seasons_m);
  fingerprint_value(hash, data.spawning_fraction_m);
  fingerprint_vector(hash, data.natural_mortality_at_age_m);
  fingerprint_vector(hash, data.maturity_at_age_m);
  fingerprint_vector(hash, data.weight_at_age_m);
  fingerprint_vector(hash, data.movement_matrix_m);
  fingerprint_vector(hash, data.regional_recruit_proportions_m);
  fingerprint_vector(hash, data.effort_m);
  fingerprint_vector(hash, data.observed_index_m);
  fingerprint_vector(hash, data.observed_retained_biomass_m);
  fingerprint_vector(hash, data.observed_discard_biomass_m);
  fingerprint_vector(hash, data.observed_total_catch_m);
  fingerprint_vector(hash, data.catch_units_m);
  fingerprint_vector(hash, data.availability_surface_m);
  fingerprint_vector(hash, data.observed_catch_numbers_m);
  fingerprint_value(hash, controls.phase_m);
  fingerprint_value(hash, controls.estimate_availability_scales_m);
  fingerprint_value(hash, controls.availability_by_fleet_only_m);
  fingerprint_value(hash, controls.share_movement_across_seasons_m);
  fingerprint_value(hash, controls.use_index_likelihood_m);
  fingerprint_value(hash, controls.use_catch_composition_likelihood_m);
  fingerprint_value(hash, controls.use_retained_biomass_likelihood_m);
  fingerprint_value(hash, controls.use_discard_biomass_likelihood_m);
  fingerprint_value(hash, controls.use_catch_conditioning_m);
  fingerprint_value(hash, controls.use_priors_m);
  fingerprint_value(hash, controls.composition_likelihood_weight_m);
  fingerprint_value(hash, controls.sd_prior_log_q_m);
  fingerprint_value(hash, controls.sd_prior_log_sel50_m);
  fingerprint_value(hash, controls.sd_prior_log_sel_slope_m);
  fingerprint_value(hash, controls.sd_prior_retention50_raw_m);
  fingerprint_value(hash, controls.sd_prior_log_retention_slope_m);
  fingerprint_value(hash, controls.sd_prior_log_availability_scale_m);
  fingerprint_value(hash, controls.sd_prior_move_logit_m);
  fingerprint_value(hash, controls.sd_prior_log_sigma_m);
  fingerprint_value(hash, controls.movement_smoothing_weight_m);
  fingerprint_value(hash, controls.availability_smoothing_weight_m);
  std::ostringstream text;
  text << std::hex << std::setw(16) << std::setfill('0') << hash;
  return text.str();
}

std::string
sampling_geometry_fingerprint(const quadra::TunaSpatialAssessmentData &data,
                              const quadra::TunaAssessmentControls &controls,
                              const std::vector<double> &parameters,
                              const std::vector<size_t> &active_positions,
                              const std::vector<double> &scales) {
  std::uint64_t hash = 1469598103934665603ULL;
  const std::string assessment = assessment_fingerprint(data, controls);
  for (const char value : assessment)
    fingerprint_value(hash, value);
  fingerprint_vector(hash, parameters);
  fingerprint_vector(hash, active_positions);
  fingerprint_vector(hash, scales);
  const double curvature_step = 1e-3;
  const double eigenvalue_floor_fraction = 1e-2;
  fingerprint_value(hash, curvature_step);
  fingerprint_value(hash, eigenvalue_floor_fraction);
  std::ostringstream text;
  text << std::hex << std::setw(16) << std::setfill('0') << hash;
  return text.str();
}

std::string file_fingerprint(const std::string &path) {
  std::ifstream input(path, std::ios::binary);
  if (!input)
    throw std::runtime_error("could not read artifact for fingerprint: " +
                             path);
  std::uint64_t hash = 1469598103934665603ULL;
  char buffer[1 << 16];
  while (input) {
    input.read(buffer, sizeof(buffer));
    const std::streamsize count = input.gcount();
    for (std::streamsize i = 0; i < count; ++i) {
      hash ^= static_cast<unsigned char>(buffer[i]);
      hash *= 1099511628211ULL;
    }
  }
  std::ostringstream text;
  text << std::hex << std::setw(16) << std::setfill('0') << hash;
  return text.str();
}

std::unordered_map<std::string, std::string>
read_flow_manifest(const std::string &artifact_path) {
  const std::string manifest_path = artifact_path + ".manifest";
  std::ifstream input(manifest_path);
  if (!input)
    throw std::runtime_error("missing transport-flow compatibility manifest: " +
                             manifest_path);
  std::unordered_map<std::string, std::string> values;
  std::string line;
  while (std::getline(input, line)) {
    const size_t separator = line.find('=');
    if (separator != std::string::npos)
      values[line.substr(0, separator)] = line.substr(separator + 1);
  }
  return values;
}

void validate_flow_manifest(const std::string &artifact_path,
                            const std::vector<std::string> &parameter_names,
                            const std::string &expected_assessment_fingerprint,
                            const std::string &expected_geometry_fingerprint) {
  const auto manifest = read_flow_manifest(artifact_path);
  const auto required = [&](const std::string &key) -> const std::string & {
    const auto found = manifest.find(key);
    if (found == manifest.end() || found->second.empty())
      throw std::runtime_error(
          "transport-flow manifest missing required field '" + key +
          "': " + artifact_path + ".manifest");
    return found->second;
  };
  std::ostringstream joined_names;
  for (size_t j = 0; j < parameter_names.size(); ++j)
    joined_names << (j == 0 ? "" : ",") << parameter_names[j];
  if (required("manifest_version") != "1")
    throw std::runtime_error("unsupported transport-flow manifest version");
  if (std::stoull(required("dimension")) != parameter_names.size())
    throw std::runtime_error(
        "transport-flow dimension does not match active model");
  if (required("parameter_names") != joined_names.str())
    throw std::runtime_error(
        "transport-flow parameter names/order do not match active model");
  if (required("assessment_fingerprint") != expected_assessment_fingerprint)
    throw std::runtime_error("transport flow was trained for a different "
                             "assessment/data fingerprint");
  if (required("geometry_fingerprint") != expected_geometry_fingerprint)
    throw std::runtime_error("transport flow is stale for the current fitted "
                             "mode or whitening geometry");
  if (required("artifact_hash_fnv1a64") != file_fingerprint(artifact_path))
    throw std::runtime_error(
        "transport-flow artifact content does not match its manifest");
  required("architecture");
  required("training_draw_sources");
  required("training_draw_hashes_fnv1a64");
  required("validation_nll");
  required("training_seed");
  const double inverse_error = std::stod(required("inverse_max_error"));
  if (!std::isfinite(inverse_error) || inverse_error > 1e-4)
    throw std::runtime_error(
        "native transport-flow inverse round-trip validation failed");
}

bool load_whitening_cache(const std::string &path,
                          const std::string &expected_fingerprint,
                          size_t dimension,
                          std::vector<std::vector<double>> &whitening) {
  std::ifstream input(path);
  std::string label, fingerprint;
  size_t cached_dimension = 0;
  if (!input || !std::getline(input, label, ',') ||
      !std::getline(input, fingerprint) || label != "fingerprint" ||
      fingerprint != expected_fingerprint || !std::getline(input, label, ',') ||
      !(input >> cached_dimension) || label != "dimension" ||
      cached_dimension != dimension)
    return false;
  input.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  whitening.assign(dimension, std::vector<double>(dimension));
  for (size_t i = 0; i < dimension; ++i) {
    std::string row;
    if (!std::getline(input, row))
      return false;
    std::istringstream values(row);
    for (size_t j = 0; j < dimension; ++j) {
      if (!(values >> whitening[i][j]) || !std::isfinite(whitening[i][j]))
        return false;
      if (j + 1 < dimension && values.get() != ',')
        return false;
    }
  }
  return true;
}

void write_whitening_cache(const std::string &path,
                           const std::string &fingerprint,
                           const std::vector<std::vector<double>> &whitening) {
  std::ostringstream csv;
  csv << std::setprecision(17) << "fingerprint," << fingerprint
      << "\ndimension," << whitening.size() << "\n";
  for (const auto &row : whitening) {
    for (size_t j = 0; j < row.size(); ++j)
      csv << (j == 0 ? "" : ",") << row[j];
    csv << "\n";
  }
  quadra::write_text_file(path, csv.str());
}

void write_fit_checkpoint(
    const std::string &out_dir,
    const quadra::AdvancedSpatialTunaAssessmentModel &model,
    const quadra::TunaFitResult &fit, const std::string &fingerprint) {
  std::ostringstream parameters;
  parameters << std::setprecision(17) << "parameter,value\n";
  const auto names = model.parameter_names();
  for (size_t i = 0; i < names.size(); ++i)
    parameters << names[i] << "," << fit.best_parameters_m[i] << "\n";
  quadra::write_text_file(out_dir + "/fit_checkpoint_parameters.csv",
                          parameters.str());
  std::ostringstream metadata;
  metadata << std::setprecision(17)
           << "fingerprint,converged,gradient_norm,nll\n"
           << fingerprint << "," << (fit.converged_m ? 1 : 0) << ","
           << fit.gradient_norm_m << "," << fit.summary_m.nll_m << "\n";
  quadra::write_text_file(out_dir + "/fit_checkpoint_metadata.csv",
                          metadata.str());
}

quadra::TunaFitResult
load_fit_checkpoint(const std::string &out_dir,
                    const quadra::TunaSpatialAssessmentData &data,
                    const quadra::TunaAssessmentControls &controls,
                    const std::string &expected_fingerprint) {
  std::ifstream metadata(out_dir + "/fit_checkpoint_metadata.csv");
  std::string header;
  std::string row;
  if (!std::getline(metadata, header) || !std::getline(metadata, row))
    throw std::runtime_error("missing or invalid fit checkpoint metadata");
  std::istringstream metadata_row(row);
  std::string fingerprint, converged, gradient, nll;
  std::getline(metadata_row, fingerprint, ',');
  std::getline(metadata_row, converged, ',');
  std::getline(metadata_row, gradient, ',');
  std::getline(metadata_row, nll, ',');
  if (fingerprint != expected_fingerprint)
    throw std::runtime_error(
        "fit checkpoint fingerprint does not match current inputs");
  if (converged != "1")
    throw std::runtime_error("fit checkpoint is not converged");

  quadra::AdvancedSpatialTunaAssessmentModel model(data, controls);
  const auto expected_names = model.parameter_names();
  std::ifstream parameters(out_dir + "/fit_checkpoint_parameters.csv");
  if (!std::getline(parameters, header))
    throw std::runtime_error("missing fit checkpoint parameters");
  std::vector<double> values;
  std::string parameter_row;
  size_t index = 0;
  while (std::getline(parameters, parameter_row)) {
    const size_t comma = parameter_row.find(',');
    if (comma == std::string::npos || index >= expected_names.size() ||
        parameter_row.substr(0, comma) != expected_names[index])
      throw std::runtime_error(
          "fit checkpoint parameter names do not match model");
    values.push_back(std::stod(parameter_row.substr(comma + 1)));
    ++index;
  }
  if (values.size() != expected_names.size())
    throw std::runtime_error("fit checkpoint parameter count mismatch");
  quadra::TunaFitResult fit;
  fit.converged_m = true;
  fit.gradient_norm_m = std::stod(gradient);
  fit.best_parameters_m = std::move(values);
  fit.message_m = "loaded validated fit checkpoint";
  fit.summary_m = quadra::evaluate_at_parameters(
      data, controls, fit.best_parameters_m, "baseline");
  return fit;
}

struct StandardizedNoncenteredTunaPosterior {
  const quadra::AdvancedSpatialTunaAssessmentModel *model_m;
  size_t sigma_index_m;
  std::vector<size_t> recruitment_indices_m;
  std::vector<size_t> active_indices_m;
  std::vector<double> noncentered_mode_m;
  std::vector<double> scales_m;
  std::vector<std::vector<double>> whitening_m;

  template <class Type> Type operator()(const std::vector<Type> &q) const {
    std::vector<Type> centered(noncentered_mode_m.begin(),
                               noncentered_mode_m.end());
    for (size_t j = 0; j < active_indices_m.size(); ++j) {
      Type standardized = Type(0.0);
      for (size_t k = 0; k < active_indices_m.size(); ++k)
        standardized += Type(whitening_m[j][k]) * q[k];
      centered[active_indices_m[j]] =
          Type(noncentered_mode_m[active_indices_m[j]]) +
          Type(scales_m[j]) * standardized;
    }
    const Type sigma = exp(centered[sigma_index_m]);
    for (const size_t index : recruitment_indices_m)
      centered[index] = sigma * centered[index];
    quadra::ModelReportContext context;
    return -model_m->evaluate(centered, context) +
           Type(static_cast<double>(recruitment_indices_m.size())) * log(sigma);
  }
};

struct FastMarginalEvaluation {
  quadra::LaplaceObjectiveResult objective_m;
  std::vector<double> gradient_m;
  bool success_m = false;
  double elapsed_ms_m = 0.0;
};

struct MarginalEvaluationProfile {
  std::atomic<long long> evaluations_m{0};
  std::atomic<long long> gradient_evaluations_m{0};
  std::atomic<long long> failures_m{0};
  std::atomic<long long> newton_iterations_m{0};
  std::atomic<long long> total_us_m{0};
  std::atomic<long long> mode_solve_us_m{0};
  std::atomic<long long> objective_us_m{0};
  std::atomic<long long> exact_gradient_us_m{0};

  void reset() {
    evaluations_m = 0;
    gradient_evaluations_m = 0;
    failures_m = 0;
    newton_iterations_m = 0;
    total_us_m = 0;
    mode_solve_us_m = 0;
    objective_us_m = 0;
    exact_gradient_us_m = 0;
  }
};

MarginalEvaluationProfile g_marginal_profile;

void record_marginal_evaluation(const FastMarginalEvaluation &evaluation) {
  ++g_marginal_profile.evaluations_m;
  if (!evaluation.success_m)
    ++g_marginal_profile.failures_m;
  g_marginal_profile.newton_iterations_m +=
      evaluation.objective_m.newton_iterations_m;
  g_marginal_profile.total_us_m +=
      static_cast<long long>(std::llround(1000.0 * evaluation.elapsed_ms_m));
  g_marginal_profile.mode_solve_us_m += static_cast<long long>(
      std::llround(1000.0 * evaluation.objective_m.mode_solve_ms_m));
  g_marginal_profile.objective_us_m += static_cast<long long>(
      std::llround(1000.0 * evaluation.objective_m.total_ms_m));
  if (!evaluation.gradient_m.empty()) {
    ++g_marginal_profile.gradient_evaluations_m;
    g_marginal_profile.exact_gradient_us_m +=
        static_cast<long long>(std::llround(
            1000.0 * std::max(0.0, evaluation.elapsed_ms_m -
                                       evaluation.objective_m.total_ms_m)));
  }
}

void write_marginal_evaluation_profile(const std::string &path) {
  const long long evaluations = g_marginal_profile.evaluations_m.load();
  const long long failures = g_marginal_profile.failures_m.load();
  const long long gradient_evaluations =
      g_marginal_profile.gradient_evaluations_m.load();
  const long long iterations = g_marginal_profile.newton_iterations_m.load();
  const double total_ms = g_marginal_profile.total_us_m.load() / 1000.0;
  const double objective_ms = g_marginal_profile.objective_us_m.load() / 1000.0;
  const double mode_ms = g_marginal_profile.mode_solve_us_m.load() / 1000.0;
  const double exact_gradient_ms =
      g_marginal_profile.exact_gradient_us_m.load() / 1000.0;
  const double divisor = std::max(1LL, evaluations);
  std::ostringstream csv;
  csv << std::setprecision(17) << "metric,value\n"
      << "evaluations," << evaluations << "\n"
      << "gradient_evaluations," << gradient_evaluations << "\n"
      << "failures," << failures << "\n"
      << "newton_iterations," << iterations << "\n"
      << "total_ms," << total_ms << "\n"
      << "objective_ms," << objective_ms << "\n"
      << "mode_solve_ms," << mode_ms << "\n"
      << "exact_gradient_ms," << exact_gradient_ms << "\n"
      << "wrapper_overhead_ms,"
      << std::max(0.0, total_ms - objective_ms - exact_gradient_ms) << "\n"
      << "mean_evaluation_ms," << total_ms / divisor << "\n"
      << "mean_objective_ms," << objective_ms / divisor << "\n"
      << "mean_mode_solve_ms," << mode_ms / divisor << "\n"
      << "mean_exact_gradient_ms,"
      << exact_gradient_ms / std::max(1LL, gradient_evaluations) << "\n"
      << "mean_newton_iterations," << static_cast<double>(iterations) / divisor
      << "\n";
  quadra::write_text_file(path, csv.str());
}

class FastMarginalLaplaceEvaluator {
  struct CombinedObjective {
    quadra::AdvancedSpatialTunaAssessmentModel *model_m;
    quadra::ParameterPartition partition_m;

    template <class Type>
    Type operator()(const std::vector<Type> &fixed_random) const {
      const size_t n_fixed = partition_m.fixed_indices_m.size();
      std::vector<Type> fixed(fixed_random.begin(),
                              fixed_random.begin() + n_fixed);
      std::vector<Type> random(fixed_random.begin() + n_fixed,
                               fixed_random.end());
      std::vector<Type> full =
          quadra::merge_parameters(fixed, random, partition_m);
      quadra::ModelReportContext context;
      model_m->initialize(context);
      return model_m->evaluate(full, context);
    }
  };

public:
  FastMarginalLaplaceEvaluator(
      quadra::AdvancedSpatialTunaAssessmentModel &model,
      std::vector<double> random_mode, quadra::ParameterPartition partition,
      const quadra::LaplaceObjectiveOptions &options)
      : model_m(&model), partition_m(std::move(partition)), options_m(options),
        last_random_mode_m(std::move(random_mode)) {
    objective_m.reset(new quadra::stats::LaplaceEvaluator<
                      quadra::AdvancedSpatialTunaAssessmentModel>(
        model, last_random_mode_m, partition_m, options_m));
  }

  FastMarginalEvaluation evaluate(const std::vector<double> &fixed,
                                  bool compute_gradient = true) {
    const auto started = std::chrono::steady_clock::now();
    FastMarginalEvaluation out;
    const size_t rebuild_interval = static_cast<size_t>(
        std::max(1, configured_int("QUADRA_TUNA_NUTS_TAPE_REBUILD_INTERVAL",
                                   "sampling.tape_rebuild_interval", 1)));
    if (evaluation_count_m > 0 && evaluation_count_m % rebuild_interval == 0) {
      objective_m.reset(new quadra::stats::LaplaceEvaluator<
                        quadra::AdvancedSpatialTunaAssessmentModel>(
          *model_m, last_random_mode_m, partition_m, options_m));
    }
    ++evaluation_count_m;
    out.objective_m = objective_m->evaluate(fixed);
    if (out.objective_m.converged_m)
      last_random_mode_m = out.objective_m.u_hat_m;
    if (!out.objective_m.converged_m || !out.objective_m.logdet_ok_m) {
      out.elapsed_ms_m = std::chrono::duration<double, std::milli>(
                             std::chrono::steady_clock::now() - started)
                             .count();
      return out;
    }
    if (!compute_gradient) {
      out.success_m = true;
      out.elapsed_ms_m = std::chrono::duration<double, std::milli>(
                             std::chrono::steady_clock::now() - started)
                             .count();
      return out;
    }

    quadra::laplace::SparseHuuFactorization factor(
        out.objective_m.hessian_random_m);
    const size_t n_fixed = fixed.size();
    const size_t n_random = out.objective_m.u_hat_m.size();
    std::vector<Eigen::VectorXd> mode_directions(n_fixed);
    for (size_t j = 0; j < n_fixed; ++j)
      mode_directions[j] =
          (-factor.solve(Eigen::VectorXd(out.objective_m.mixed_hessian_m.col(
               static_cast<Eigen::Index>(j)))))
              .eval();

    Eigen::MatrixXd inverse(n_random, n_random);
    for (size_t j = 0; j < n_random; ++j) {
      Eigen::VectorXd unit = Eigen::VectorXd::Zero(n_random);
      unit[static_cast<Eigen::Index>(j)] = 1.0;
      inverse.col(static_cast<Eigen::Index>(j)) = factor.solve(unit);
    }
    inverse = 0.5 * (inverse + inverse.transpose()).eval();
    Eigen::LLT<Eigen::MatrixXd> inverse_factor(inverse);
    if (inverse_factor.info() != Eigen::Success) {
      out.elapsed_ms_m = std::chrono::duration<double, std::milli>(
                             std::chrono::steady_clock::now() - started)
                             .count();
      return out;
    }
    const Eigen::MatrixXd root = inverse_factor.matrixL();

    std::vector<double> combined = fixed;
    combined.insert(combined.end(), out.objective_m.u_hat_m.begin(),
                    out.objective_m.u_hat_m.end());
    out.gradient_m = out.objective_m.gradient_fixed_joint_m;
    std::vector<double> trace_terms(n_fixed, 0.0);
    const size_t worker_count = std::min<size_t>(
        n_fixed, static_cast<size_t>(std::max(
                     1, configured_int("QUADRA_TUNA_GRADIENT_WORKERS",
                                       "sampling.gradient_workers", 4))));
    const auto compute_trace_partition = [&](size_t worker) {
      // Each concurrent worker owns a model copy because model
      // initialization maintains evaluation-local state.
      auto local_model = *model_m;
      CombinedObjective local_objective{&local_model, partition_m};
      for (size_t j = worker; j < n_fixed; j += worker_count) {
        std::vector<double> total(combined.size(), 0.0);
        total[j] = 1.0;
        for (size_t r = 0; r < n_random; ++r)
          total[n_fixed + r] = mode_directions[j][static_cast<Eigen::Index>(r)];
        const double base =
            had::third_directional_derivative(local_objective, combined, total);
        double trace = 0.0;
        for (size_t column = 0; column < n_random; ++column) {
          std::vector<double> plus = total;
          std::vector<double> minus = total;
          for (size_t r = 0; r < n_random; ++r) {
            const double twice_root =
                2.0 * root(static_cast<Eigen::Index>(r),
                           static_cast<Eigen::Index>(column));
            plus[n_fixed + r] += twice_root;
            minus[n_fixed + r] -= twice_root;
          }
          trace += (had::third_directional_derivative(local_objective, combined,
                                                      plus) +
                    had::third_directional_derivative(local_objective, combined,
                                                      minus) -
                    2.0 * base) /
                   24.0;
        }
        trace_terms[j] = trace;
      }
    };
    if (worker_count == 1) {
      // Avoid creating and destroying one nested thread for every
      // leapfrog evaluation in the normal multi-chain workflow.
      compute_trace_partition(0);
    } else {
      std::vector<std::future<void>> workers;
      workers.reserve(worker_count);
      for (size_t worker = 0; worker < worker_count; ++worker)
        workers.emplace_back(
            std::async(std::launch::async, compute_trace_partition, worker));
      for (auto &worker : workers)
        worker.get();
    }
    for (size_t j = 0; j < n_fixed; ++j)
      out.gradient_m[j] += 0.5 * trace_terms[j];
    out.success_m = true;
    out.elapsed_ms_m = std::chrono::duration<double, std::milli>(
                           std::chrono::steady_clock::now() - started)
                           .count();
    return out;
  }

private:
  quadra::AdvancedSpatialTunaAssessmentModel *model_m;
  quadra::ParameterPartition partition_m;
  quadra::LaplaceObjectiveOptions options_m;
  std::vector<double> last_random_mode_m;
  size_t evaluation_count_m = 0;
  std::unique_ptr<quadra::stats::LaplaceEvaluator<
      quadra::AdvancedSpatialTunaAssessmentModel>>
      objective_m;
};

class MarginalLaplaceTunaPosterior {
public:
  MarginalLaplaceTunaPosterior(const quadra::TunaSpatialAssessmentData &data,
                               const quadra::TunaAssessmentControls &controls,
                               const quadra::ParameterPartition &partition,
                               std::vector<double> fixed_mode,
                               std::vector<double> random_mode,
                               std::vector<size_t> active_fixed_positions,
                               std::vector<double> scales,
                               std::vector<std::vector<double>> whitening,
                               bool require_gradient = true)
      : model_m(new quadra::AdvancedSpatialTunaAssessmentModel(data, controls)),
        fixed_mode_m(std::move(fixed_mode)),
        random_mode_m(std::move(random_mode)),
        active_fixed_positions_m(std::move(active_fixed_positions)),
        scales_m(std::move(scales)), whitening_m(std::move(whitening)) {
    quadra::LaplaceObjectiveOptions objective_options;
    objective_options.include_constant_m = true;
    objective_options.compute_mixed_derivatives_m = require_gradient;
    objective_options.newton_m.max_iterations_m = 300;
    objective_options.newton_m.gradient_tolerance_m = 1e-5;
    objective_options.newton_m.step_tolerance_m = 1e-10;
    objective_options.newton_m.initial_step_scale_m =
        configured_double("QUADRA_TUNA_NUTS_NEWTON_INITIAL_STEP",
                          "sampling.newton_initial_step", 1.0);
    objective_options.newton_m.min_step_scale_m = 1e-10;
    objective_options.newton_m.sufficient_decrease_m = 1e-6;
    objective_options.newton_m.hessian_drop_tol_m = 0.0;
    objective_options.newton_m.use_backtracking_m = true;
    evaluator_m.reset(new FastMarginalLaplaceEvaluator(
        *model_m, random_mode_m, partition, objective_options));
  }

  MarginalLaplaceTunaPosterior(MarginalLaplaceTunaPosterior &&) = default;
  MarginalLaplaceTunaPosterior &
  operator=(MarginalLaplaceTunaPosterior &&) = default;
  MarginalLaplaceTunaPosterior(const MarginalLaplaceTunaPosterior &) = delete;
  MarginalLaplaceTunaPosterior &
  operator=(const MarginalLaplaceTunaPosterior &) = delete;

  template <class Type> Type operator()(const std::vector<Type> &q) const {
    std::vector<double> fixed = fixed_mode_m;
    for (size_t j = 0; j < q.size(); ++j) {
      double transformed = 0.0;
      for (size_t k = 0; k < q.size(); ++k)
        transformed += whitening_m[j][k] * quadra::value_of(q[k]);
      fixed[active_fixed_positions_m[j]] += scales_m[j] * transformed;
    }

    // ExactLaplaceEvaluator uses its own AD tapes. Restore the outer
    // sampling tape before attaching the supplied marginal gradient.
    auto *outer_graph = had::g_ADGraph;
    constexpr bool type_requires_gradient =
        !std::is_same<typename std::decay<Type>::type, double>::value;
    const FastMarginalEvaluation evaluated =
        evaluator_m->evaluate(fixed, type_requires_gradient);
    had::g_ADGraph = outer_graph;
    record_marginal_evaluation(evaluated);
    if (!evaluated.success_m)
      return Type(-std::numeric_limits<double>::infinity());
    random_mode_m = evaluated.objective_m.u_hat_m;

    Type log_density(-evaluated.objective_m.laplace_objective_m);
    if constexpr (type_requires_gradient) {
      for (size_t j = 0; j < q.size(); ++j) {
        const double q_value = quadra::value_of(q[j]);
        double gradient = 0.0;
        for (size_t i = 0; i < q.size(); ++i)
          gradient -= evaluated.gradient_m[active_fixed_positions_m[i]] *
                      scales_m[i] * whitening_m[i][j];
        log_density += Type(gradient) * (q[j] - Type(q_value));
      }
    }
    return log_density;
  }

private:
  std::unique_ptr<quadra::AdvancedSpatialTunaAssessmentModel> model_m;
  std::vector<double> fixed_mode_m;
  mutable std::vector<double> random_mode_m;
  std::vector<size_t> active_fixed_positions_m;
  std::vector<double> scales_m;
  std::vector<std::vector<double>> whitening_m;
  mutable std::unique_ptr<FastMarginalLaplaceEvaluator> evaluator_m;
};

double squared_norm(const std::vector<double> &x) {
  double out = 0.0;
  for (const double value : x)
    out += value * value;
  return out;
}

template <class Target>
quadra::sampling::NutsResult
sample_pcn_chain(Target &target, const std::vector<double> &initial, int warmup,
                 int samples, double initial_beta, double target_acceptance,
                 bool adapt_beta, std::uint64_t seed) {
  if (!(initial_beta > 0.0 && initial_beta <= 1.0) ||
      !(target_acceptance > 0.0 && target_acceptance < 1.0))
    throw std::invalid_argument("invalid pCN sampler configuration");
  std::mt19937_64 rng(seed);
  std::normal_distribution<double> normal;
  std::uniform_real_distribution<double> uniform;
  std::vector<double> current = initial;
  double current_log_density = target(current);
  if (!std::isfinite(current_log_density))
    throw std::runtime_error("pCN initial state has non-finite density");
  double current_residual = current_log_density + 0.5 * squared_norm(current);
  double beta = initial_beta;
  double logit_beta = std::log(beta / (1.0 - std::min(beta, 1.0 - 1e-12)));
  double warmup_acceptance = 0.0;
  double sampling_acceptance = 0.0;
  quadra::sampling::NutsResult out;
  out.draws.reserve(static_cast<size_t>(samples));
  out.log_density.reserve(static_cast<size_t>(samples));
  out.energy.reserve(static_cast<size_t>(samples));
  for (int iteration = 0; iteration < warmup + samples; ++iteration) {
    const double persistence = std::sqrt(std::max(0.0, 1.0 - beta * beta));
    std::vector<double> proposal(current.size());
    for (size_t j = 0; j < current.size(); ++j)
      proposal[j] = persistence * current[j] + beta * normal(rng);
    const double proposed_log_density = target(proposal);
    double acceptance = 0.0;
    double proposed_residual = -std::numeric_limits<double>::infinity();
    if (std::isfinite(proposed_log_density)) {
      proposed_residual = proposed_log_density + 0.5 * squared_norm(proposal);
      acceptance = std::min(
          1.0, std::exp(std::min(0.0, proposed_residual - current_residual)));
    }
    if (uniform(rng) < acceptance) {
      current = std::move(proposal);
      current_log_density = proposed_log_density;
      current_residual = proposed_residual;
    }
    if (iteration < warmup) {
      warmup_acceptance += acceptance;
      if (adapt_beta) {
        const double gain = 0.5 / std::sqrt(static_cast<double>(iteration + 1));
        logit_beta += gain * (acceptance - target_acceptance);
        beta = 1.0 / (1.0 + std::exp(-logit_beta));
        beta = std::min(0.999, std::max(0.01, beta));
      }
    } else {
      sampling_acceptance += acceptance;
      out.draws.push_back(current);
      out.log_density.push_back(current_log_density);
      out.energy.push_back(-current_log_density);
    }
  }
  out.diagnostics.step_size = beta;
  out.diagnostics.warmup_mean_acceptance = warmup_acceptance / warmup;
  out.diagnostics.mean_acceptance = sampling_acceptance / samples;
  double energy_mean = 0.0;
  for (double value : out.energy)
    energy_mean += value;
  energy_mean /= out.energy.size();
  double variance = 0.0;
  double difference = 0.0;
  for (size_t i = 0; i < out.energy.size(); ++i) {
    variance += (out.energy[i] - energy_mean) * (out.energy[i] - energy_mean);
    if (i > 0)
      difference += (out.energy[i] - out.energy[i - 1]) *
                    (out.energy[i] - out.energy[i - 1]);
  }
  out.diagnostics.energy_bfmi = variance > 0.0 ? difference / variance : 1.0;
  return out;
}

template <class TargetFactory>
quadra::sampling::NutsWorkflowResult
run_pcn_workflow(TargetFactory target_factory, const std::vector<double> &mode,
                 std::vector<std::string> names,
                 const quadra::sampling::NutsWorkflowOptions &options,
                 double beta, double target_acceptance, bool adapt_beta) {
  quadra::sampling::validate_parameter_names(names, mode.size());
  quadra::sampling::NutsWorkflowResult out;
  out.parameter_names = std::move(names);
  out.options = options;
  out.initial_states = quadra::sampling::initialize_nuts_chains(
      mode, options.chains, options.initial_jitter,
      options.initialization_seed);
  out.fit.chains.resize(options.chains);
  std::vector<std::future<quadra::sampling::NutsResult>> futures;
  for (size_t chain = 0; chain < options.chains; ++chain)
    futures.push_back(std::async(std::launch::async, [&, chain] {
      auto target = target_factory(chain);
      return sample_pcn_chain(
          target, out.initial_states[chain], options.sampler.warmup,
          options.sampler.samples, beta, target_acceptance, adapt_beta,
          options.sampler.seed + UINT64_C(0x9e3779b97f4a7c15) * (chain + 1));
    }));
  for (size_t chain = 0; chain < options.chains; ++chain)
    out.fit.chains[chain] = futures[chain].get();
  out.fit.diagnostics =
      quadra::sampling::compute_multi_chain_diagnostics(out.fit.chains);
  out.health = quadra::sampling::assess_nuts_health(out.fit, options.health);
  return out;
}

double log_add_exp(double left, double right);

struct StudentTProposal {
  Eigen::VectorXd mean_m;
  Eigen::MatrixXd root_m;
  Eigen::MatrixXd precision_m;
  double df_m = 8.0;
  double log_normalizer_m = 0.0;

  double log_density(const std::vector<double> &value) const {
    const Eigen::Map<const Eigen::VectorXd> mapped(
        value.data(), static_cast<Eigen::Index>(value.size()));
    const Eigen::VectorXd difference = mapped - mean_m;
    const double quadratic = difference.dot(precision_m * difference);
    return log_normalizer_m -
           0.5 * (df_m + value.size()) * std::log1p(quadratic / df_m);
  }

  std::vector<double> draw(std::mt19937_64 &rng) const {
    std::normal_distribution<double> normal;
    std::gamma_distribution<double> gamma(df_m / 2.0, 2.0);
    Eigen::VectorXd z(mean_m.size());
    for (Eigen::Index j = 0; j < z.size(); ++j)
      z[j] = normal(rng);
    const double scale = std::sqrt(df_m / gamma(rng));
    const Eigen::VectorXd sampled = mean_m + scale * root_m * z;
    return std::vector<double>(sampled.data(), sampled.data() + sampled.size());
  }
};

struct GaussianMixtureProposal {
  struct Component {
    Eigen::VectorXd mean_m;
    Eigen::MatrixXd root_m;
    Eigen::MatrixXd precision_m;
    double weight_m = 0.0;
    double log_normalizer_m = 0.0;
  };
  std::vector<Component> components_m;

  double log_density(const std::vector<double> &value) const {
    const Eigen::Map<const Eigen::VectorXd> mapped(
        value.data(), static_cast<Eigen::Index>(value.size()));
    double total = -std::numeric_limits<double>::infinity();
    for (const auto &component : components_m) {
      const Eigen::VectorXd difference = mapped - component.mean_m;
      const double term =
          std::log(component.weight_m) + component.log_normalizer_m -
          0.5 * difference.dot(component.precision_m * difference);
      total = std::isfinite(total) ? log_add_exp(total, term) : term;
    }
    return total;
  }

  std::vector<double> draw(std::mt19937_64 &rng) const {
    std::uniform_real_distribution<double> uniform;
    std::normal_distribution<double> normal;
    double selected = uniform(rng);
    const Component *component = &components_m.back();
    for (const auto &candidate : components_m) {
      selected -= candidate.weight_m;
      if (selected <= 0.0) {
        component = &candidate;
        break;
      }
    }
    Eigen::VectorXd z(component->mean_m.size());
    for (Eigen::Index j = 0; j < z.size(); ++j)
      z[j] = normal(rng);
    const Eigen::VectorXd sampled = component->mean_m + component->root_m * z;
    return std::vector<double>(sampled.data(), sampled.data() + sampled.size());
  }
};

struct GaussianKdeProposal {
  std::vector<Eigen::VectorXd> centers_m;
  Eigen::MatrixXd root_m;
  Eigen::MatrixXd precision_m;
  double log_component_normalizer_m = 0.0;

  double log_density(const std::vector<double> &value) const {
    const Eigen::Map<const Eigen::VectorXd> mapped(
        value.data(), static_cast<Eigen::Index>(value.size()));
    double total = -std::numeric_limits<double>::infinity();
    for (const auto &center : centers_m) {
      const Eigen::VectorXd difference = mapped - center;
      const double term = log_component_normalizer_m -
                          0.5 * difference.dot(precision_m * difference);
      total = std::isfinite(total) ? log_add_exp(total, term) : term;
    }
    return total - std::log(static_cast<double>(centers_m.size()));
  }

  std::vector<double> draw(std::mt19937_64 &rng) const {
    std::uniform_int_distribution<size_t> select(0, centers_m.size() - 1);
    std::normal_distribution<double> normal;
    Eigen::VectorXd z(root_m.cols());
    for (Eigen::Index j = 0; j < z.size(); ++j)
      z[j] = normal(rng);
    const Eigen::VectorXd sampled = centers_m[select(rng)] + root_m * z;
    return std::vector<double>(sampled.data(), sampled.data() + sampled.size());
  }
};

struct PolynomialFlowProposal {
  Eigen::VectorXd mean_m;
  Eigen::MatrixXd root_m;
  Eigen::MatrixXd inverse_root_m;
  std::vector<Eigen::VectorXd> coefficients_m;
  Eigen::VectorXd sigma_m;
  double log_linear_jacobian_m = 0.0;
  double training_nll_m = 0.0;
  double validation_nll_m = 0.0;

  static Eigen::VectorXd features(const Eigen::VectorXd &x, Eigen::Index j) {
    const Eigen::Index paired = std::min<Eigen::Index>(j, 6);
    const Eigen::Index interactions = paired * (paired - 1) / 2;
    Eigen::VectorXd out(1 + 2 * j + interactions);
    Eigen::Index position = 0;
    out[position++] = 1.0;
    for (Eigen::Index i = 0; i < j; ++i)
      out[position++] = x[i];
    for (Eigen::Index i = 0; i < j; ++i)
      out[position++] = x[i] * x[i];
    for (Eigen::Index left = 0; left < paired; ++left)
      for (Eigen::Index right = left + 1; right < paired; ++right)
        out[position++] = x[left] * x[right];
    return out;
  }

  double log_density(const std::vector<double> &value) const {
    const Eigen::Map<const Eigen::VectorXd> q(
        value.data(), static_cast<Eigen::Index>(value.size()));
    const Eigen::VectorXd x = inverse_root_m * (q - mean_m);
    double out = log_linear_jacobian_m;
    for (Eigen::Index j = 0; j < x.size(); ++j) {
      const double conditional_mean =
          coefficients_m[static_cast<size_t>(j)].dot(features(x, j));
      const double residual = (x[j] - conditional_mean) / sigma_m[j];
      out += -0.5 * residual * residual - std::log(sigma_m[j]) -
             0.5 * std::log(2.0 * M_PI);
    }
    return out;
  }

  std::vector<double> draw(std::mt19937_64 &rng) const {
    std::normal_distribution<double> normal;
    Eigen::VectorXd x = Eigen::VectorXd::Zero(mean_m.size());
    for (Eigen::Index j = 0; j < x.size(); ++j)
      x[j] = coefficients_m[static_cast<size_t>(j)].dot(features(x, j)) +
             sigma_m[j] * normal(rng);
    const Eigen::VectorXd q = mean_m + root_m * x;
    return std::vector<double>(q.data(), q.data() + q.size());
  }
};

struct QFlowProposal {
  quadra::transport::DependencyFreeRealNVP native_m;
  size_t dimension_m = 0;
  Eigen::VectorXd fixed_origin_m;
  Eigen::MatrixXd q_to_fixed_m;
  Eigen::MatrixXd fixed_to_q_m;

  QFlowProposal(const std::string &path, const std::vector<double> &fixed_mode,
                const std::vector<size_t> &active_positions,
                const std::vector<double> &scales,
                const std::vector<std::vector<double>> &whitening,
                const std::vector<std::string> &parameter_names,
                const std::string &assessment_fingerprint,
                const std::string &geometry_fingerprint)
      : dimension_m(active_positions.size()),
        fixed_origin_m(active_positions.size()),
        q_to_fixed_m(active_positions.size(), active_positions.size()) {
    validate_flow_manifest(path, parameter_names, assessment_fingerprint,
                           geometry_fingerprint);
    const auto manifest = read_flow_manifest(path);
    const std::string architecture = manifest.at("architecture");
    if (architecture.rfind("qflow_realnvp_cpp;", 0) != 0)
      throw std::runtime_error(
          "transport artifact is not a dependency-free QFLOW archive");
    native_m.load(path);
    if (native_m.dimension() != dimension_m)
      throw std::runtime_error("QFLOW archive dimension mismatch");
    for (size_t i = 0; i < dimension_m; ++i) {
      fixed_origin_m[static_cast<Eigen::Index>(i)] =
          fixed_mode[active_positions[i]];
      for (size_t j = 0; j < dimension_m; ++j)
        q_to_fixed_m(static_cast<Eigen::Index>(i),
                     static_cast<Eigen::Index>(j)) =
            scales[i] * whitening[i][j];
    }
    fixed_to_q_m = q_to_fixed_m.inverse();
  }

  double log_density(const std::vector<double> &value) const {
    const Eigen::Map<const Eigen::VectorXd> q(
        value.data(), static_cast<Eigen::Index>(dimension_m));
    const Eigen::VectorXd fixed = fixed_origin_m + q_to_fixed_m * q;
    return native_m.log_density(
        std::vector<double>(fixed.data(), fixed.data() + fixed.size()));
  }

  std::vector<double> draw(std::mt19937_64 &rng) const {
    std::normal_distribution<double> normal;
    std::vector<double> noise(dimension_m);
    for (double &value : noise)
      value = normal(rng);
    const std::vector<double> values = native_m.inverse(noise);
    Eigen::VectorXd fixed(dimension_m);
    for (size_t j = 0; j < dimension_m; ++j)
      fixed[static_cast<Eigen::Index>(j)] = values[j];
    const Eigen::VectorXd q = fixed_to_q_m * (fixed - fixed_origin_m);
    return std::vector<double>(q.data(), q.data() + q.size());
  }
};

struct QFlowMixtureProposal {
  std::vector<QFlowProposal> components_m;

  QFlowMixtureProposal(const std::vector<std::string> &paths,
                       const std::vector<double> &fixed_mode,
                       const std::vector<size_t> &active_positions,
                       const std::vector<double> &scales,
                       const std::vector<std::vector<double>> &whitening,
                       const std::vector<std::string> &parameter_names,
                       const std::string &assessment_fingerprint,
                       const std::string &geometry_fingerprint) {
    if (paths.empty())
      throw std::invalid_argument("transport flow mixture is empty");
    components_m.reserve(paths.size());
    for (const auto &path : paths)
      components_m.emplace_back(path, fixed_mode, active_positions, scales,
                                whitening, parameter_names,
                                assessment_fingerprint, geometry_fingerprint);
  }

  double log_density(const std::vector<double> &value) const {
    double total = -std::numeric_limits<double>::infinity();
    for (const auto &component : components_m) {
      const double term = component.log_density(value);
      total = std::isfinite(total) ? log_add_exp(total, term) : term;
    }
    return total - std::log(static_cast<double>(components_m.size()));
  }

  std::vector<double> draw(std::mt19937_64 &rng) const {
    std::uniform_int_distribution<size_t> select(0, components_m.size() - 1);
    return components_m[select(rng)].draw(rng);
  }
};

std::vector<std::vector<double>>
read_training_draws(const std::string &path,
                    const std::vector<std::string> &names,
                    const std::vector<double> &fixed_mode,
                    const std::vector<size_t> &active_positions,
                    const std::vector<double> &scales,
                    const std::vector<std::vector<double>> &whitening) {
  std::ifstream input(path);
  std::string line;
  if (!input || !std::getline(input, line))
    throw std::runtime_error("could not read proposal training draws: " + path);
  std::unordered_map<std::string, size_t> position;
  for (size_t j = 0; j < names.size(); ++j)
    position[names[j]] = j;
  std::vector<std::vector<double>> fixed_draws;
  std::vector<double> current(names.size());
  std::vector<bool> seen(names.size(), false);
  std::string previous_key;
  auto finish = [&]() {
    if (previous_key.empty())
      return;
    if (std::find(seen.begin(), seen.end(), false) != seen.end())
      throw std::runtime_error("proposal training draw is incomplete");
    fixed_draws.push_back(current);
    std::fill(seen.begin(), seen.end(), false);
  };
  while (std::getline(input, line)) {
    std::istringstream row(line);
    std::string chain, iteration, parameter, value;
    std::getline(row, chain, ',');
    std::getline(row, iteration, ',');
    std::getline(row, parameter, ',');
    std::getline(row, value, ',');
    const std::string key = chain + ":" + iteration;
    if (!previous_key.empty() && key != previous_key)
      finish();
    previous_key = key;
    const auto found = position.find(parameter);
    if (found == position.end())
      throw std::runtime_error("proposal training parameter mismatch: " +
                               parameter);
    current[found->second] = std::stod(value);
    seen[found->second] = true;
  }
  finish();
  if (fixed_draws.size() < 2 * names.size())
    throw std::runtime_error("too few proposal training draws");

  const size_t dimension = names.size();
  Eigen::MatrixXd transform(dimension, dimension);
  for (size_t i = 0; i < dimension; ++i)
    for (size_t j = 0; j < dimension; ++j)
      transform(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) =
          whitening[i][j];
  Eigen::FullPivLU<Eigen::MatrixXd> inverse(transform);
  if (!inverse.isInvertible())
    throw std::runtime_error("whitening transform is not invertible");
  std::vector<std::vector<double>> standardized;
  standardized.reserve(fixed_draws.size());
  for (const auto &draw : fixed_draws) {
    Eigen::VectorXd displaced(dimension);
    for (size_t j = 0; j < dimension; ++j)
      displaced[static_cast<Eigen::Index>(j)] =
          (draw[j] - fixed_mode[active_positions[j]]) / scales[j];
    const Eigen::VectorXd q = inverse.solve(displaced);
    standardized.emplace_back(q.data(), q.data() + q.size());
  }
  return standardized;
}

StudentTProposal
fit_student_t_proposal(const std::vector<std::vector<double>> &draws, double df,
                       double inflation) {
  if (!(df > 2.0) || !(inflation > 0.0))
    throw std::invalid_argument("invalid Student-t proposal configuration");
  const Eigen::Index dimension =
      static_cast<Eigen::Index>(draws.front().size());
  Eigen::VectorXd mean = Eigen::VectorXd::Zero(dimension);
  for (const auto &draw : draws)
    mean += Eigen::Map<const Eigen::VectorXd>(draw.data(), dimension);
  mean /= static_cast<double>(draws.size());
  Eigen::MatrixXd covariance = Eigen::MatrixXd::Zero(dimension, dimension);
  for (const auto &draw : draws) {
    const Eigen::VectorXd difference =
        Eigen::Map<const Eigen::VectorXd>(draw.data(), dimension) - mean;
    covariance.noalias() += difference * difference.transpose();
  }
  covariance /= static_cast<double>(draws.size() - 1);
  covariance *= inflation * (df - 2.0) / df;
  covariance.diagonal().array() +=
      std::max(1e-10, covariance.diagonal().mean() * 1e-8);
  Eigen::LLT<Eigen::MatrixXd> factor(covariance);
  if (factor.info() != Eigen::Success)
    throw std::runtime_error(
        "Student-t proposal covariance is not positive definite");
  StudentTProposal out;
  out.mean_m = mean;
  out.root_m = factor.matrixL();
  out.precision_m =
      factor.solve(Eigen::MatrixXd::Identity(dimension, dimension));
  out.df_m = df;
  double log_determinant = 0.0;
  for (Eigen::Index j = 0; j < dimension; ++j)
    log_determinant += 2.0 * std::log(out.root_m(j, j));
  out.log_normalizer_m =
      std::lgamma(0.5 * (df + dimension)) - std::lgamma(0.5 * df) -
      0.5 * (dimension * std::log(df * M_PI) + log_determinant);
  return out;
}

GaussianMixtureProposal
fit_gaussian_mixture_proposal(const std::vector<std::vector<double>> &draws,
                              size_t component_count, double inflation,
                              double shrinkage) {
  if (component_count < 2 || component_count > draws.size() / 10 ||
      !(inflation > 0.0) || !(shrinkage >= 0.0 && shrinkage < 1.0))
    throw std::invalid_argument(
        "invalid Gaussian-mixture transport configuration");
  const Eigen::Index dimension =
      static_cast<Eigen::Index>(draws.front().size());
  std::vector<Eigen::VectorXd> points;
  points.reserve(draws.size());
  for (const auto &draw : draws)
    points.emplace_back(
        Eigen::Map<const Eigen::VectorXd>(draw.data(), dimension));
  Eigen::VectorXd global_mean = Eigen::VectorXd::Zero(dimension);
  for (const auto &point : points)
    global_mean += point;
  global_mean /= static_cast<double>(points.size());
  Eigen::MatrixXd global_covariance =
      Eigen::MatrixXd::Zero(dimension, dimension);
  for (const auto &point : points) {
    const Eigen::VectorXd difference = point - global_mean;
    global_covariance.noalias() += difference * difference.transpose();
  }
  global_covariance /= static_cast<double>(points.size() - 1);

  // Deterministic farthest-point initialization followed by Lloyd steps.
  std::vector<Eigen::VectorXd> centers{points.front()};
  while (centers.size() < component_count) {
    size_t farthest = 0;
    double farthest_distance = -1.0;
    for (size_t i = 0; i < points.size(); ++i) {
      double nearest = std::numeric_limits<double>::infinity();
      for (const auto &center : centers)
        nearest = std::min(nearest, (points[i] - center).squaredNorm());
      if (nearest > farthest_distance) {
        farthest_distance = nearest;
        farthest = i;
      }
    }
    centers.push_back(points[farthest]);
  }
  std::vector<size_t> assignment(points.size());
  for (int iteration = 0; iteration < 30; ++iteration) {
    bool changed = false;
    std::vector<Eigen::VectorXd> sums(component_count,
                                      Eigen::VectorXd::Zero(dimension));
    std::vector<size_t> counts(component_count, 0);
    for (size_t i = 0; i < points.size(); ++i) {
      size_t nearest_component = 0;
      double nearest_distance = std::numeric_limits<double>::infinity();
      for (size_t component = 0; component < component_count; ++component) {
        const double distance = (points[i] - centers[component]).squaredNorm();
        if (distance < nearest_distance) {
          nearest_distance = distance;
          nearest_component = component;
        }
      }
      changed = changed || iteration == 0 || assignment[i] != nearest_component;
      assignment[i] = nearest_component;
      sums[nearest_component] += points[i];
      ++counts[nearest_component];
    }
    for (size_t component = 0; component < component_count; ++component)
      if (counts[component] > 0)
        centers[component] =
            sums[component] / static_cast<double>(counts[component]);
    if (!changed)
      break;
  }

  GaussianMixtureProposal out;
  for (size_t component = 0; component < component_count; ++component) {
    std::vector<size_t> members;
    for (size_t i = 0; i < assignment.size(); ++i)
      if (assignment[i] == component)
        members.push_back(i);
    if (members.size() <= static_cast<size_t>(dimension))
      continue;
    Eigen::VectorXd mean = Eigen::VectorXd::Zero(dimension);
    for (size_t index : members)
      mean += points[index];
    mean /= static_cast<double>(members.size());
    Eigen::MatrixXd covariance = Eigen::MatrixXd::Zero(dimension, dimension);
    for (size_t index : members) {
      const Eigen::VectorXd difference = points[index] - mean;
      covariance.noalias() += difference * difference.transpose();
    }
    covariance /= static_cast<double>(members.size() - 1);
    covariance = inflation * ((1.0 - shrinkage) * covariance +
                              shrinkage * global_covariance);
    covariance.diagonal().array() +=
        std::max(1e-10, covariance.diagonal().mean() * 1e-8);
    Eigen::LLT<Eigen::MatrixXd> factor(covariance);
    if (factor.info() != Eigen::Success)
      throw std::runtime_error("Gaussian-mixture component covariance failed");
    GaussianMixtureProposal::Component fitted;
    fitted.mean_m = mean;
    fitted.root_m = factor.matrixL();
    fitted.precision_m =
        factor.solve(Eigen::MatrixXd::Identity(dimension, dimension));
    fitted.weight_m = static_cast<double>(members.size()) /
                      static_cast<double>(points.size());
    double log_determinant = 0.0;
    for (Eigen::Index j = 0; j < dimension; ++j)
      log_determinant += 2.0 * std::log(fitted.root_m(j, j));
    fitted.log_normalizer_m =
        -0.5 * (dimension * std::log(2.0 * M_PI) + log_determinant);
    out.components_m.push_back(std::move(fitted));
  }
  double total_weight = 0.0;
  for (const auto &component : out.components_m)
    total_weight += component.weight_m;
  if (!(total_weight > 0.0))
    throw std::runtime_error("Gaussian-mixture transport has no components");
  for (auto &component : out.components_m)
    component.weight_m /= total_weight;
  return out;
}

GaussianKdeProposal
fit_gaussian_kde_proposal(const std::vector<std::vector<double>> &draws,
                          double bandwidth) {
  if (!(bandwidth > 0.0))
    throw std::invalid_argument("invalid KDE transport bandwidth");
  const StudentTProposal covariance =
      fit_student_t_proposal(draws, 100.0, 100.0 / 98.0);
  GaussianKdeProposal out;
  out.root_m = bandwidth * covariance.root_m;
  const Eigen::Index dimension = out.root_m.rows();
  Eigen::LLT<Eigen::MatrixXd> factor(out.root_m * out.root_m.transpose());
  out.precision_m =
      factor.solve(Eigen::MatrixXd::Identity(dimension, dimension));
  double log_determinant = 0.0;
  for (Eigen::Index j = 0; j < dimension; ++j)
    log_determinant += 2.0 * std::log(out.root_m(j, j));
  out.log_component_normalizer_m =
      -0.5 * (dimension * std::log(2.0 * M_PI) + log_determinant);
  out.centers_m.reserve(draws.size());
  for (const auto &draw : draws)
    out.centers_m.emplace_back(
        Eigen::Map<const Eigen::VectorXd>(draw.data(), dimension));
  return out;
}

PolynomialFlowProposal
fit_polynomial_flow_proposal(const std::vector<std::vector<double>> &draws,
                             double ridge) {
  if (!(ridge >= 0.0))
    throw std::invalid_argument("invalid polynomial-flow ridge penalty");
  const Eigen::Index dimension =
      static_cast<Eigen::Index>(draws.front().size());
  PolynomialFlowProposal out;
  out.mean_m = Eigen::VectorXd::Zero(dimension);
  for (const auto &draw : draws)
    out.mean_m += Eigen::Map<const Eigen::VectorXd>(draw.data(), dimension);
  out.mean_m /= static_cast<double>(draws.size());
  Eigen::MatrixXd covariance = Eigen::MatrixXd::Zero(dimension, dimension);
  for (const auto &draw : draws) {
    const Eigen::VectorXd difference =
        Eigen::Map<const Eigen::VectorXd>(draw.data(), dimension) - out.mean_m;
    covariance.noalias() += difference * difference.transpose();
  }
  covariance /= static_cast<double>(draws.size() - 1);
  covariance.diagonal().array() +=
      std::max(1e-10, covariance.diagonal().mean() * 1e-6);
  Eigen::LLT<Eigen::MatrixXd> factor(covariance);
  if (factor.info() != Eigen::Success)
    throw std::runtime_error("polynomial-flow whitening failed");
  out.root_m = factor.matrixL();
  out.inverse_root_m = out.root_m.inverse();
  out.log_linear_jacobian_m = 0.0;
  for (Eigen::Index j = 0; j < dimension; ++j)
    out.log_linear_jacobian_m -= std::log(out.root_m(j, j));
  std::vector<Eigen::VectorXd> x;
  x.reserve(draws.size());
  for (const auto &draw : draws)
    x.push_back(out.inverse_root_m *
                (Eigen::Map<const Eigen::VectorXd>(draw.data(), dimension) -
                 out.mean_m));
  out.coefficients_m.resize(static_cast<size_t>(dimension));
  out.sigma_m.resize(dimension);
  for (Eigen::Index j = 0; j < dimension; ++j) {
    const Eigen::Index feature_count =
        PolynomialFlowProposal::features(x.front(), j).size();
    Eigen::MatrixXd gram = Eigen::MatrixXd::Zero(feature_count, feature_count);
    Eigen::VectorXd cross = Eigen::VectorXd::Zero(feature_count);
    size_t training_count = 0;
    for (size_t i = 0; i < x.size(); ++i)
      if (i % 5 != 0) {
        const Eigen::VectorXd feature =
            PolynomialFlowProposal::features(x[i], j);
        gram.noalias() += feature * feature.transpose();
        cross.noalias() += feature * x[i][j];
        ++training_count;
      }
    for (Eigen::Index k = 1; k < feature_count; ++k)
      gram(k, k) += ridge * static_cast<double>(training_count);
    Eigen::LDLT<Eigen::MatrixXd> solver(gram);
    if (solver.info() != Eigen::Success)
      throw std::runtime_error("polynomial-flow regression failed");
    out.coefficients_m[static_cast<size_t>(j)] = solver.solve(cross);
    double residual_sum = 0.0;
    for (size_t i = 0; i < x.size(); ++i)
      if (i % 5 != 0) {
        const double residual =
            x[i][j] - out.coefficients_m[static_cast<size_t>(j)].dot(
                          PolynomialFlowProposal::features(x[i], j));
        residual_sum += residual * residual;
      }
    out.sigma_m[j] = std::max(
        0.15, std::sqrt(residual_sum / static_cast<double>(training_count)));
  }
  size_t training_count = 0;
  size_t validation_count = 0;
  for (size_t i = 0; i < draws.size(); ++i)
    if (i % 5 == 0) {
      out.validation_nll_m -= out.log_density(draws[i]);
      ++validation_count;
    } else {
      out.training_nll_m -= out.log_density(draws[i]);
      ++training_count;
    }
  out.training_nll_m /= static_cast<double>(training_count);
  out.validation_nll_m /= static_cast<double>(validation_count);
  return out;
}

double log_add_exp(double left, double right) {
  const double maximum = std::max(left, right);
  return maximum +
         std::log(std::exp(left - maximum) + std::exp(right - maximum));
}

double local_normal_log_density(const std::vector<double> &to,
                                const std::vector<double> &from, double scale) {
  double quadratic = 0.0;
  for (size_t j = 0; j < to.size(); ++j)
    quadratic += (to[j] - from[j]) * (to[j] - from[j]);
  return -0.5 * quadratic / (scale * scale) -
         static_cast<double>(to.size()) *
             (std::log(scale) + 0.5 * std::log(2.0 * M_PI));
}

template <class Target, class Proposal>
quadra::sampling::NutsResult
sample_independence_chain(Target &target, const std::vector<double> &initial,
                          int warmup, int samples, const Proposal &proposal,
                          double global_weight, double local_scale,
                          std::uint64_t seed) {
  if (!(global_weight > 0.0 && global_weight <= 1.0) || !(local_scale > 0.0))
    throw std::invalid_argument("invalid independence sampler mixture");
  std::mt19937_64 rng(seed);
  std::normal_distribution<double> normal;
  std::uniform_real_distribution<double> uniform;
  std::vector<double> current = initial;
  double current_target = target(current);
  if (!std::isfinite(current_target))
    throw std::runtime_error(
        "independence sampler initial density is not finite");
  double acceptance_sum = 0.0;
  double warmup_acceptance = 0.0;
  quadra::sampling::NutsResult out;
  out.draws.reserve(samples);
  out.log_density.reserve(samples);
  out.energy.reserve(samples);
  const auto mixture_log_density = [&](const std::vector<double> &to,
                                       const std::vector<double> &from) {
    const double global = std::log(global_weight) + proposal.log_density(to);
    if (global_weight >= 1.0)
      return global;
    const double local = std::log1p(-global_weight) +
                         local_normal_log_density(to, from, local_scale);
    return log_add_exp(global, local);
  };
  for (int iteration = 0; iteration < warmup + samples; ++iteration) {
    std::vector<double> proposed;
    if (uniform(rng) < global_weight)
      proposed = proposal.draw(rng);
    else {
      proposed = current;
      for (double &value : proposed)
        value += local_scale * normal(rng);
    }
    const double proposed_target = target(proposed);
    double acceptance = 0.0;
    if (std::isfinite(proposed_target))
      acceptance = std::min(
          1.0,
          std::exp(std::min(0.0, proposed_target - current_target +
                                     mixture_log_density(current, proposed) -
                                     mixture_log_density(proposed, current))));
    if (uniform(rng) < acceptance) {
      current = std::move(proposed);
      current_target = proposed_target;
    }
    if (iteration < warmup)
      warmup_acceptance += acceptance;
    else {
      acceptance_sum += acceptance;
      out.draws.push_back(current);
      out.log_density.push_back(current_target);
      out.energy.push_back(-current_target);
    }
  }
  out.diagnostics.mean_acceptance = acceptance_sum / samples;
  out.diagnostics.warmup_mean_acceptance = warmup_acceptance / warmup;
  out.diagnostics.step_size = local_scale;
  out.diagnostics.energy_bfmi = 1.0;
  return out;
}

template <class TargetFactory, class Proposal>
quadra::sampling::NutsWorkflowResult run_independence_workflow(
    TargetFactory target_factory, const std::vector<double> &mode,
    std::vector<std::string> names,
    const quadra::sampling::NutsWorkflowOptions &options,
    const Proposal &proposal, double global_weight, double local_scale,
    const std::vector<std::vector<double>> *training_draws = nullptr) {
  quadra::sampling::NutsWorkflowResult out;
  out.parameter_names = std::move(names);
  out.options = options;
  out.initial_states = quadra::sampling::initialize_nuts_chains(
      mode, options.chains, options.initial_jitter,
      options.initialization_seed);
  if (training_draws != nullptr)
    for (size_t chain = 0; chain < options.chains; ++chain)
      out.initial_states[chain] =
          (*training_draws)[(chain + 1) * training_draws->size() /
                            (options.chains + 1)];
  out.fit.chains.resize(options.chains);
  std::vector<std::future<quadra::sampling::NutsResult>> futures;
  for (size_t chain = 0; chain < options.chains; ++chain)
    futures.push_back(std::async(std::launch::async, [&, chain] {
      auto target = target_factory(chain);
      return sample_independence_chain(
          target, out.initial_states[chain], options.sampler.warmup,
          options.sampler.samples, proposal, global_weight, local_scale,
          options.sampler.seed + UINT64_C(0x9e3779b97f4a7c15) * (chain + 1));
    }));
  for (size_t chain = 0; chain < options.chains; ++chain)
    out.fit.chains[chain] = futures[chain].get();
  out.fit.diagnostics =
      quadra::sampling::compute_multi_chain_diagnostics(out.fit.chains);
  out.health = quadra::sampling::assess_nuts_health(out.fit, options.health);
  return out;
}

template <class Target, class Proposal>
quadra::sampling::NutsResult
sample_isir_chain(Target &target, const std::vector<double> &initial,
                  int warmup, int samples, const Proposal &proposal,
                  size_t candidates, std::uint64_t seed) {
  if (candidates < 2)
    throw std::invalid_argument("i-SIR requires at least two candidates");
  std::mt19937_64 rng(seed);
  std::uniform_real_distribution<double> uniform;
  std::vector<double> current = initial;
  double current_target = target(current);
  double move_sum = 0.0;
  double warmup_move_sum = 0.0;
  quadra::sampling::NutsResult out;
  for (int iteration = 0; iteration < warmup + samples; ++iteration) {
    std::vector<std::vector<double>> states(candidates);
    std::vector<double> targets(candidates);
    std::vector<double> log_weights(candidates);
    states[0] = current;
    targets[0] = current_target;
    log_weights[0] = current_target - proposal.log_density(current);
    for (size_t candidate = 1; candidate < candidates; ++candidate) {
      states[candidate] = proposal.draw(rng);
      targets[candidate] = target(states[candidate]);
      log_weights[candidate] =
          std::isfinite(targets[candidate])
              ? targets[candidate] - proposal.log_density(states[candidate])
              : -std::numeric_limits<double>::infinity();
    }
    const double maximum =
        *std::max_element(log_weights.begin(), log_weights.end());
    double total = 0.0;
    for (double &weight : log_weights) {
      weight = std::exp(weight - maximum);
      total += weight;
    }
    double selected_weight = uniform(rng) * total;
    size_t selected = candidates - 1;
    for (size_t candidate = 0; candidate < candidates; ++candidate) {
      selected_weight -= log_weights[candidate];
      if (selected_weight <= 0.0) {
        selected = candidate;
        break;
      }
    }
    const double moved = selected == 0 ? 0.0 : 1.0;
    current = std::move(states[selected]);
    current_target = targets[selected];
    if (iteration < warmup)
      warmup_move_sum += moved;
    else {
      move_sum += moved;
      out.draws.push_back(current);
      out.log_density.push_back(current_target);
      out.energy.push_back(-current_target);
    }
  }
  out.diagnostics.mean_acceptance = move_sum / samples;
  out.diagnostics.warmup_mean_acceptance = warmup_move_sum / warmup;
  out.diagnostics.step_size = static_cast<double>(candidates);
  out.diagnostics.energy_bfmi = 1.0;
  return out;
}

template <class TargetFactory, class Proposal>
quadra::sampling::NutsWorkflowResult
run_isir_workflow(TargetFactory target_factory, std::vector<std::string> names,
                  const quadra::sampling::NutsWorkflowOptions &options,
                  const Proposal &proposal, size_t candidates,
                  const std::vector<std::vector<double>> &training_draws) {
  quadra::sampling::NutsWorkflowResult out;
  out.parameter_names = std::move(names);
  out.options = options;
  out.initial_states.resize(options.chains);
  for (size_t chain = 0; chain < options.chains; ++chain)
    out.initial_states[chain] =
        training_draws[(chain + 1) * training_draws.size() /
                       (options.chains + 1)];
  out.fit.chains.resize(options.chains);
  std::vector<std::future<quadra::sampling::NutsResult>> futures;
  for (size_t chain = 0; chain < options.chains; ++chain)
    futures.push_back(std::async(std::launch::async, [&, chain] {
      auto target = target_factory(chain);
      return sample_isir_chain(
          target, out.initial_states[chain], options.sampler.warmup,
          options.sampler.samples, proposal, candidates,
          options.sampler.seed + UINT64_C(0x9e3779b97f4a7c15) * (chain + 1));
    }));
  for (size_t chain = 0; chain < options.chains; ++chain)
    out.fit.chains[chain] = futures[chain].get();
  out.fit.diagnostics =
      quadra::sampling::compute_multi_chain_diagnostics(out.fit.chains);
  out.health = quadra::sampling::assess_nuts_health(out.fit, options.health);
  return out;
}

template <class Target>
quadra::sampling::NutsResult
sample_covariance_rw_chain(Target &target, const std::vector<double> &initial,
                           int warmup, int samples, const Eigen::MatrixXd &root,
                           double scale, std::uint64_t seed) {
  std::mt19937_64 rng(seed);
  std::normal_distribution<double> normal;
  std::uniform_real_distribution<double> uniform;
  std::vector<double> current = initial;
  double current_target = target(current);
  double warmup_acceptance = 0.0;
  double sampling_acceptance = 0.0;
  quadra::sampling::NutsResult out;
  for (int iteration = 0; iteration < warmup + samples; ++iteration) {
    Eigen::VectorXd z(root.cols());
    for (Eigen::Index j = 0; j < z.size(); ++j)
      z[j] = normal(rng);
    const Eigen::VectorXd increment = scale * root * z;
    std::vector<double> proposed = current;
    for (size_t j = 0; j < proposed.size(); ++j)
      proposed[j] += increment[static_cast<Eigen::Index>(j)];
    const double proposed_target = target(proposed);
    const double acceptance =
        std::isfinite(proposed_target)
            ? std::min(1.0, std::exp(std::min(0.0, proposed_target -
                                                       current_target)))
            : 0.0;
    if (uniform(rng) < acceptance) {
      current = std::move(proposed);
      current_target = proposed_target;
    }
    if (iteration < warmup)
      warmup_acceptance += acceptance;
    else {
      sampling_acceptance += acceptance;
      out.draws.push_back(current);
      out.log_density.push_back(current_target);
      out.energy.push_back(-current_target);
    }
  }
  out.diagnostics.mean_acceptance = sampling_acceptance / samples;
  out.diagnostics.warmup_mean_acceptance = warmup_acceptance / warmup;
  out.diagnostics.step_size = scale;
  out.diagnostics.energy_bfmi = 1.0;
  return out;
}

template <class TargetFactory>
quadra::sampling::NutsWorkflowResult run_covariance_rw_workflow(
    TargetFactory target_factory, const std::vector<double> &,
    std::vector<std::string> names,
    const quadra::sampling::NutsWorkflowOptions &options,
    const Eigen::MatrixXd &root, double scale,
    const std::vector<std::vector<double>> &training_draws) {
  quadra::sampling::NutsWorkflowResult out;
  out.parameter_names = std::move(names);
  out.options = options;
  out.initial_states.resize(options.chains);
  for (size_t chain = 0; chain < options.chains; ++chain)
    out.initial_states[chain] =
        training_draws[(chain + 1) * training_draws.size() /
                       (options.chains + 1)];
  out.fit.chains.resize(options.chains);
  std::vector<std::future<quadra::sampling::NutsResult>> futures;
  for (size_t chain = 0; chain < options.chains; ++chain)
    futures.push_back(std::async(std::launch::async, [&, chain] {
      auto target = target_factory(chain);
      return sample_covariance_rw_chain(
          target, out.initial_states[chain], options.sampler.warmup,
          options.sampler.samples, root, scale,
          options.sampler.seed + UINT64_C(0x9e3779b97f4a7c15) * (chain + 1));
    }));
  for (size_t chain = 0; chain < options.chains; ++chain)
    out.fit.chains[chain] = futures[chain].get();
  out.fit.diagnostics =
      quadra::sampling::compute_multi_chain_diagnostics(out.fit.chains);
  out.health = quadra::sampling::assess_nuts_health(out.fit, options.health);
  return out;
}

void run_optional_marginal_nuts(const std::string &out_dir,
                                const quadra::TunaSpatialAssessmentData &data,
                                const quadra::TunaAssessmentControls &controls,
                                const quadra::TunaFitResult &fit) {
  if (!configured_bool("QUADRA_TUNA_RUN_NUTS", "sampling.enabled", false))
    return;
  if (!fit.converged_m)
    throw std::runtime_error(
        "marginal AD-NUTS requires a converged assessment mode");

  quadra::AdvancedSpatialTunaAssessmentModel prototype(data, controls);
  const quadra::ParameterSet parameters = prototype.parameter_set();
  const quadra::ParameterPartition partition =
      quadra::partition_parameters(parameters);
  const auto split = quadra::split_parameters(fit.best_parameters_m, partition);
  const std::vector<std::string> all_names = prototype.parameter_names();
  const int anchor_fleet =
      configured_int("QUADRA_TUNA_ANCHOR_FLEET", "model.anchor_fleet", 1);
  const std::vector<std::string> locked_names = {
      "logit_steepness", "log_index_q_fleet_" + std::to_string(anchor_fleet)};
  std::vector<size_t> active_positions;
  std::vector<std::string> names;
  std::vector<double> scales;
  for (size_t j = 0; j < partition.fixed_indices_m.size(); ++j) {
    const std::string &name = all_names[partition.fixed_indices_m[j]];
    if (std::find(locked_names.begin(), locked_names.end(), name) !=
        locked_names.end())
      continue;
    active_positions.push_back(j);
    names.push_back(name);
    double scale = 0.25;
    if (name.find("sel50_raw_") == 0 || name.find("retention50_raw_") == 0)
      scale = 0.5;
    else if (name.find("move_logit_") == 0)
      scale = 0.4;
    scales.push_back(scale);
  }

  const size_t dimension = active_positions.size();
  const std::string geometry_fingerprint = sampling_geometry_fingerprint(
      data, controls, fit.best_parameters_m, active_positions, scales);
  const std::string geometry_cache_path =
      out_dir + "/marginal_whitening_cache.csv";
  std::vector<std::vector<double>> whitening;
  const bool geometry_cache_hit = load_whitening_cache(
      geometry_cache_path, geometry_fingerprint, dimension, whitening);

  // Validate the reduced exact trace identity against Quadra's general
  // dense third-order evaluator at the fitted mode. This audit runs once
  // per workflow, not once per transition.
  if (!geometry_cache_hit) {
    quadra::LaplaceObjectiveOptions audit_objective_options;
    audit_objective_options.include_constant_m = true;
    audit_objective_options.compute_mixed_derivatives_m = true;
    audit_objective_options.newton_m.max_iterations_m = 300;
    audit_objective_options.newton_m.gradient_tolerance_m = 1e-5;
    audit_objective_options.newton_m.step_tolerance_m = 1e-10;
    audit_objective_options.newton_m.initial_step_scale_m = 0.20;
    audit_objective_options.newton_m.min_step_scale_m = 1e-10;
    audit_objective_options.newton_m.sufficient_decrease_m = 1e-6;
    audit_objective_options.newton_m.hessian_drop_tol_m = 0.0;
    audit_objective_options.newton_m.use_backtracking_m = true;
    quadra::AdvancedSpatialTunaAssessmentModel audit_model(data, controls);
    FastMarginalLaplaceEvaluator fast_audit(audit_model, split.random_m,
                                            partition, audit_objective_options);
    const FastMarginalEvaluation fast_result =
        fast_audit.evaluate(split.fixed_m);
    quadra::laplace::ExactLaplaceGradientEngineOptions audit_engine_options;
    audit_engine_options.discover_active_directions = false;
    audit_engine_options.stream_dense_hdot_trace = true;
    audit_engine_options.hdot_workers = 4;
    quadra::stats::ExactLaplaceEvaluator<
        quadra::AdvancedSpatialTunaAssessmentModel>
        reference_audit(audit_model, split.fixed_m, split.random_m, partition,
                        audit_objective_options, audit_engine_options);
    const quadra::stats::ExactLaplaceResult reference_result =
        reference_audit.evaluate(split.fixed_m, split.random_m);
    if (!fast_result.success_m || !reference_result.success)
      throw std::runtime_error("marginal gradient audit evaluation failed");
    double maximum_gradient_difference = 0.0;
    std::ostringstream gradient_audit;
    gradient_audit << std::setprecision(17)
                   << "parameter,reference_gradient,reduced_trace_gradient,"
                      "absolute_difference\n";
    for (size_t j = 0; j < split.fixed_m.size(); ++j) {
      const double difference =
          std::abs(reference_result.gradient[j] - fast_result.gradient_m[j]);
      maximum_gradient_difference =
          std::max(maximum_gradient_difference, difference);
      gradient_audit << all_names[partition.fixed_indices_m[j]] << ","
                     << reference_result.gradient[j] << ","
                     << fast_result.gradient_m[j] << "," << difference << "\n";
    }
    quadra::write_text_file(out_dir + "/marginal_gradient_audit.csv",
                            gradient_audit.str());
    std::ostringstream performance_audit;
    performance_audit << std::setprecision(17) << "implementation,elapsed_ms\n"
                      << "general_dense_trace,"
                      << reference_result.timings.total_ms << "\n"
                      << "factored_inverse_trace," << fast_result.elapsed_ms_m
                      << "\n"
                      << "speedup,"
                      << reference_result.timings.total_ms /
                             fast_result.elapsed_ms_m
                      << "\n";
    quadra::write_text_file(out_dir + "/marginal_gradient_performance.csv",
                            performance_audit.str());
    if (maximum_gradient_difference > 1e-7)
      throw std::runtime_error(
          "reduced-trace marginal gradient failed exact audit");

    // Whiten the active fixed effects using the local Hessian of the
    // exact marginal posterior. This captures the strong midpoint/slope,
    // movement, and catchability rotations before NUTS warmup.
    Eigen::MatrixXd curvature(dimension, dimension);
    constexpr double curvature_step = 1e-3;
    for (size_t j = 0; j < dimension; ++j) {
      std::vector<double> plus = split.fixed_m;
      std::vector<double> minus = split.fixed_m;
      plus[active_positions[j]] += scales[j] * curvature_step;
      minus[active_positions[j]] -= scales[j] * curvature_step;
      const FastMarginalEvaluation plus_result = fast_audit.evaluate(plus);
      const FastMarginalEvaluation minus_result = fast_audit.evaluate(minus);
      if (!plus_result.success_m || !minus_result.success_m)
        throw std::runtime_error(
            "marginal Hessian whitening evaluation failed");
      for (size_t i = 0; i < dimension; ++i)
        curvature(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) =
            scales[i] *
            (plus_result.gradient_m[active_positions[i]] -
             minus_result.gradient_m[active_positions[i]]) /
            (2.0 * curvature_step);
    }
    curvature = 0.5 * (curvature + curvature.transpose()).eval();
    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(curvature);
    if (eigen.info() != Eigen::Success)
      throw std::runtime_error(
          "marginal Hessian whitening eigendecomposition failed");
    const double maximum_eigenvalue =
        std::max(1.0, eigen.eigenvalues().maxCoeff());
    const double eigenvalue_floor = maximum_eigenvalue * 1e-2;
    Eigen::VectorXd regularized = eigen.eigenvalues();
    size_t regularized_directions = 0;
    for (Eigen::Index j = 0; j < regularized.size(); ++j)
      if (regularized[j] < eigenvalue_floor) {
        regularized[j] = eigenvalue_floor;
        ++regularized_directions;
      }
    const Eigen::MatrixXd whitening_matrix =
        eigen.eigenvectors() *
        regularized.array().rsqrt().matrix().asDiagonal() *
        eigen.eigenvectors().transpose();
    whitening.assign(dimension, std::vector<double>(dimension));
    for (size_t i = 0; i < dimension; ++i)
      for (size_t j = 0; j < dimension; ++j)
        whitening[i][j] = whitening_matrix(static_cast<Eigen::Index>(i),
                                           static_cast<Eigen::Index>(j));
    std::ostringstream whitening_audit;
    whitening_audit << std::setprecision(17)
                    << "direction,raw_eigenvalue,regularized_eigenvalue\n";
    for (Eigen::Index j = 0; j < regularized.size(); ++j)
      whitening_audit << j + 1 << "," << eigen.eigenvalues()[j] << ","
                      << regularized[j] << "\n";
    whitening_audit << "regularized_direction_count,," << regularized_directions
                    << "\n";
    quadra::write_text_file(out_dir + "/marginal_hessian_whitening.csv",
                            whitening_audit.str());
    write_whitening_cache(geometry_cache_path, geometry_fingerprint, whitening);
  }

  // Exclude the one-time exact-gradient audit from sampler profiling.
  g_marginal_profile.reset();

  quadra::sampling::NutsWorkflowOptions options;
  options.chains = static_cast<size_t>(std::max(
      2, configured_int("QUADRA_TUNA_NUTS_CHAINS", "sampling.chains", 4)));
  options.sampler.warmup = std::max(
      1, configured_int("QUADRA_TUNA_NUTS_WARMUP", "sampling.warmup", 500));
  options.sampler.samples = std::max(
      1, configured_int("QUADRA_TUNA_NUTS_SAMPLES", "sampling.samples", 500));
  options.sampler.max_tree_depth =
      std::max(1, configured_int("QUADRA_TUNA_NUTS_MAX_TREE_DEPTH",
                                 "sampling.max_tree_depth", 10));
  options.sampler.target_acceptance = configured_double(
      "QUADRA_TUNA_NUTS_TARGET_ACCEPTANCE", "sampling.target_acceptance", 0.8);
  options.sampler.initial_step_size = configured_double(
      "QUADRA_TUNA_NUTS_INITIAL_STEP_SIZE", "sampling.initial_step_size", 0.0);
  options.sampler.adapt_diagonal_mass = configured_bool(
      "QUADRA_TUNA_NUTS_ADAPT_MASS", "sampling.adapt_mass", true);
  options.sampler.adapt_dense_mass =
      options.sampler.adapt_diagonal_mass &&
      configured_bool("QUADRA_TUNA_NUTS_DENSE_METRIC", "sampling.dense_metric",
                      true);
  options.sampler.reuse_ad_graph = false;
  options.parallel = configured_bool("QUADRA_TUNA_NUTS_PARALLEL_CHAINS",
                                     "sampling.parallel_chains", true);
  options.sampler.seed = 20260820;
  options.initialization_seed = 20260820;
  options.health.max_rhat = 1.05;

  std::vector<double> zero(active_positions.size(), 0.0);
  std::string sampler_method =
      configured_string("QUADRA_TUNA_SAMPLER", "sampling.method", "pcn");
  std::transform(sampler_method.begin(), sampler_method.end(),
                 sampler_method.begin(),
                 [](unsigned char c) { return std::tolower(c); });
  // BFMI diagnoses Hamiltonian energy transitions and has no validity
  // criterion for Metropolis pCN. R-hat and ESS remain enforced.
  if (sampler_method == "pcn" || sampler_method == "independence_t" ||
      sampler_method == "covariance_rw" || sampler_method == "transport_gmm" ||
      sampler_method == "transport_kde" || sampler_method == "transport_poly" ||
      sampler_method == "transport_flow" || sampler_method == "transport_isir")
    options.health.min_bfmi = 0.0;
  const bool sampler_requires_gradient =
      sampler_method != "pcn" && sampler_method != "independence_t" &&
      sampler_method != "covariance_rw" && sampler_method != "transport_gmm" &&
      sampler_method != "transport_kde" && sampler_method != "transport_poly" &&
      sampler_method != "transport_flow" && sampler_method != "transport_isir";
  auto target_factory = [&data, &controls, partition, fixed = split.fixed_m,
                         random = split.random_m, active_positions, scales,
                         whitening, sampler_requires_gradient](size_t) {
    return MarginalLaplaceTunaPosterior(data, controls, partition, fixed,
                                        random, active_positions, scales,
                                        whitening, sampler_requires_gradient);
  };
  quadra::sampling::NutsWorkflowResult workflow;
  if (sampler_method == "pcn")
    workflow = run_pcn_workflow(
        target_factory, zero, names, options,
        configured_double("QUADRA_TUNA_PCN_BETA", "sampling.pcn_beta", 0.5),
        configured_double("QUADRA_TUNA_PCN_TARGET_ACCEPTANCE",
                          "sampling.pcn_target_acceptance", 0.3),
        configured_bool("QUADRA_TUNA_PCN_ADAPT_BETA", "sampling.pcn_adapt_beta",
                        false));
  else if (sampler_method == "independence_t") {
    const auto training_draws = read_training_draws(
        configured_string(
            "QUADRA_TUNA_PROPOSAL_DRAWS", "sampling.proposal_draws",
            "build/assessment_outputs/benchmark_data/posterior_draws.csv"),
        names, split.fixed_m, active_positions, scales, whitening);
    const StudentTProposal proposal = fit_student_t_proposal(
        training_draws,
        configured_double("QUADRA_TUNA_PROPOSAL_DF", "sampling.proposal_df",
                          8.0),
        configured_double("QUADRA_TUNA_PROPOSAL_INFLATION",
                          "sampling.proposal_inflation", 1.2));
    workflow = run_independence_workflow(
        target_factory, zero, names, options, proposal,
        configured_double("QUADRA_TUNA_PROPOSAL_GLOBAL_WEIGHT",
                          "sampling.proposal_global_weight", 0.9),
        configured_double("QUADRA_TUNA_PROPOSAL_LOCAL_SCALE",
                          "sampling.proposal_local_scale", 0.2),
        &training_draws);
  } else if (sampler_method == "transport_gmm") {
    const auto training_draws = read_training_draws(
        configured_string(
            "QUADRA_TUNA_PROPOSAL_DRAWS", "sampling.proposal_draws",
            "build/assessment_outputs/benchmark_data/posterior_draws.csv"),
        names, split.fixed_m, active_positions, scales, whitening);
    const GaussianMixtureProposal proposal = fit_gaussian_mixture_proposal(
        training_draws,
        static_cast<size_t>(
            std::max(2, configured_int("QUADRA_TUNA_TRANSPORT_COMPONENTS",
                                       "sampling.transport_components", 6))),
        configured_double("QUADRA_TUNA_TRANSPORT_INFLATION",
                          "sampling.transport_inflation", 1.1),
        configured_double("QUADRA_TUNA_TRANSPORT_SHRINKAGE",
                          "sampling.transport_shrinkage", 0.15));
    workflow = run_independence_workflow(
        target_factory, zero, names, options, proposal,
        configured_double("QUADRA_TUNA_PROPOSAL_GLOBAL_WEIGHT",
                          "sampling.proposal_global_weight", 0.9),
        configured_double("QUADRA_TUNA_PROPOSAL_LOCAL_SCALE",
                          "sampling.proposal_local_scale", 0.2),
        &training_draws);
  } else if (sampler_method == "transport_kde") {
    const auto training_draws = read_training_draws(
        configured_string(
            "QUADRA_TUNA_PROPOSAL_DRAWS", "sampling.proposal_draws",
            "build/assessment_outputs/benchmark_data/posterior_draws.csv"),
        names, split.fixed_m, active_positions, scales, whitening);
    const GaussianKdeProposal proposal = fit_gaussian_kde_proposal(
        training_draws,
        configured_double("QUADRA_TUNA_TRANSPORT_BANDWIDTH",
                          "sampling.transport_bandwidth", 0.15));
    workflow = run_independence_workflow(
        target_factory, zero, names, options, proposal,
        configured_double("QUADRA_TUNA_PROPOSAL_GLOBAL_WEIGHT",
                          "sampling.proposal_global_weight", 0.9),
        configured_double("QUADRA_TUNA_PROPOSAL_LOCAL_SCALE",
                          "sampling.proposal_local_scale", 0.2),
        &training_draws);
  } else if (sampler_method == "transport_poly") {
    const auto training_draws = read_training_draws(
        configured_string(
            "QUADRA_TUNA_PROPOSAL_DRAWS", "sampling.proposal_draws",
            "build/assessment_outputs/benchmark_data/posterior_draws.csv"),
        names, split.fixed_m, active_positions, scales, whitening);
    const PolynomialFlowProposal proposal = fit_polynomial_flow_proposal(
        training_draws, configured_double("QUADRA_TUNA_TRANSPORT_RIDGE",
                                          "sampling.transport_ridge", 0.001));
    std::ostringstream validation;
    validation << std::setprecision(17) << "metric,value\ntraining_mean_nll,"
               << proposal.training_nll_m << "\nvalidation_mean_nll,"
               << proposal.validation_nll_m << "\ngeneralization_gap,"
               << proposal.validation_nll_m - proposal.training_nll_m << "\n";
    quadra::write_text_file(out_dir + "/transport_validation.csv",
                            validation.str());
    workflow = run_independence_workflow(
        target_factory, zero, names, options, proposal, 1.0,
        configured_double("QUADRA_TUNA_PROPOSAL_LOCAL_SCALE",
                          "sampling.proposal_local_scale", 0.2),
        &training_draws);
  } else if (sampler_method == "transport_flow" ||
             sampler_method == "transport_isir") {
    const auto training_draws = read_training_draws(
        configured_string(
            "QUADRA_TUNA_PROPOSAL_DRAWS", "sampling.proposal_draws",
            "build/assessment_outputs/benchmark_data/posterior_draws.csv"),
        names, split.fixed_m, active_positions, scales, whitening);
    const std::string model_list = configured_string(
        "QUADRA_TUNA_TRANSPORT_MODELS", "sampling.transport_models",
        configured_string(
            "QUADRA_TUNA_TRANSPORT_MODEL", "sampling.transport_model",
            "build/transport_flow/native_realnvp_20260823.qflow"));
    std::vector<std::string> model_paths;
    std::istringstream models(model_list);
    std::string model_path;
    while (std::getline(models, model_path, ','))
      if (!model_path.empty())
        model_paths.push_back(model_path);
    const QFlowMixtureProposal proposal(
        model_paths, split.fixed_m, active_positions, scales, whitening, names,
        assessment_fingerprint(data, controls), geometry_fingerprint);
    if (sampler_method == "transport_isir")
      workflow = run_isir_workflow(
          target_factory, names, options, proposal,
          static_cast<size_t>(
              std::max(2, configured_int("QUADRA_TUNA_TRANSPORT_CANDIDATES",
                                         "sampling.transport_candidates", 4))),
          training_draws);
    else
      workflow = run_independence_workflow(
          target_factory, zero, names, options, proposal,
          configured_double("QUADRA_TUNA_PROPOSAL_GLOBAL_WEIGHT",
                            "sampling.proposal_global_weight", 0.9),
          configured_double("QUADRA_TUNA_PROPOSAL_LOCAL_SCALE",
                            "sampling.proposal_local_scale", 0.2),
          &training_draws);
  } else if (sampler_method == "covariance_rw") {
    const auto training_draws = read_training_draws(
        configured_string(
            "QUADRA_TUNA_PROPOSAL_DRAWS", "sampling.proposal_draws",
            "build/assessment_outputs/benchmark_data/posterior_draws.csv"),
        names, split.fixed_m, active_positions, scales, whitening);
    const StudentTProposal covariance =
        fit_student_t_proposal(training_draws, 100.0, 100.0 / 98.0);
    workflow = run_covariance_rw_workflow(
        target_factory, zero, names, options, covariance.root_m,
        configured_double("QUADRA_TUNA_PROPOSAL_RW_SCALE",
                          "sampling.proposal_rw_scale", 0.5),
        training_draws);
  } else if (sampler_method == "nuts" || sampler_method == "ad-nuts")
    workflow = quadra::sampling::run_nuts_workflow(target_factory, zero, names,
                                                   options);
  else
    throw std::invalid_argument(
        "sampling.method must be transport_isir, transport_flow, "
        "transport_poly, transport_kde, transport_gmm, covariance_rw, "
        "independence_t, pcn, or nuts, got: " +
        sampler_method);

  write_marginal_evaluation_profile(out_dir + "/marginal_sampler_profile.csv");

  for (auto &chain : workflow.fit.chains)
    for (auto &draw : chain.draws) {
      const std::vector<double> standardized = draw;
      for (size_t j = 0; j < draw.size(); ++j) {
        double transformed = 0.0;
        for (size_t k = 0; k < draw.size(); ++k)
          transformed += whitening[j][k] * standardized[k];
        draw[j] = split.fixed_m[active_positions[j]] + scales[j] * transformed;
      }
    }
  quadra::sampling::write_posterior_draws_csv(out_dir + "/posterior_draws.csv",
                                              workflow);
  quadra::sampling::write_parameter_diagnostics_csv(
      out_dir + "/posterior_parameter_diagnostics.csv", workflow);
  quadra::sampling::write_chain_diagnostics_csv(
      out_dir + "/posterior_chain_diagnostics.csv", workflow);
  quadra::sampling::write_nuts_summary_csv(out_dir + "/sampler_summary.csv",
                                           workflow);
  quadra::sampling::write_nuts_summary_csv(out_dir + "/nuts_summary.csv",
                                           workflow);
  const int evaluations_per_transition =
      sampler_method == "transport_isir"
          ? std::max(1, configured_int("QUADRA_TUNA_TRANSPORT_CANDIDATES",
                                       "sampling.transport_candidates", 4) -
                            1)
          : 1;
  std::ostringstream sampler_identity;
  sampler_identity << "metric,value\nmethod," << sampler_method
                   << "\nevaluations_per_transition,"
                   << evaluations_per_transition << "\n"
                   << "geometry_cache_hit," << (geometry_cache_hit ? 1 : 0)
                   << "\n";
  quadra::write_text_file(out_dir + "/sampler_identity.csv",
                          sampler_identity.str());
  struct ReconstructedDraw {
    bool valid_m = false;
    std::vector<double> mode_m;
    std::vector<double> random_m;
    std::vector<double> full_m;
  };
  std::vector<std::vector<ReconstructedDraw>> reconstructed(
      workflow.fit.chains.size());
  std::vector<std::future<std::vector<ReconstructedDraw>>> futures;
  futures.reserve(workflow.fit.chains.size());
  for (size_t chain = 0; chain < workflow.fit.chains.size(); ++chain) {
    futures.emplace_back(std::async(std::launch::async, [&, chain] {
      quadra::AdvancedSpatialTunaAssessmentModel reconstruction_model(data,
                                                                      controls);
      quadra::LaplaceObjectiveOptions reconstruction_options;
      reconstruction_options.include_constant_m = true;
      reconstruction_options.compute_mixed_derivatives_m = false;
      reconstruction_options.newton_m.max_iterations_m = 300;
      reconstruction_options.newton_m.gradient_tolerance_m = 1e-5;
      reconstruction_options.newton_m.step_tolerance_m = 1e-10;
      reconstruction_options.newton_m.initial_step_scale_m =
          configured_double("QUADRA_TUNA_NUTS_NEWTON_INITIAL_STEP",
                            "sampling.newton_initial_step", 1.0);
      reconstruction_options.newton_m.min_step_scale_m = 1e-10;
      reconstruction_options.newton_m.sufficient_decrease_m = 1e-6;
      reconstruction_options.newton_m.hessian_drop_tol_m = 0.0;
      reconstruction_options.newton_m.use_backtracking_m = true;
      std::vector<double> random_initial = split.random_m;
      std::mt19937_64 rng(options.sampler.seed + 9000001u + 100003u * chain);
      std::normal_distribution<double> normal(0.0, 1.0);
      std::vector<ReconstructedDraw> result;
      result.reserve(workflow.fit.chains[chain].draws.size());
      for (const auto &draw : workflow.fit.chains[chain].draws) {
        ReconstructedDraw item;
        std::vector<double> fixed = split.fixed_m;
        for (size_t j = 0; j < active_positions.size(); ++j)
          fixed[active_positions[j]] = draw[j];
        try {
          const auto objective = quadra::evaluate_laplace_objective(
              reconstruction_model, fixed, random_initial, partition,
              reconstruction_options);
          if (!objective.converged_m || !objective.logdet_ok_m ||
              !std::isfinite(objective.laplace_objective_m)) {
            result.push_back(std::move(item));
            continue;
          }
          random_initial = objective.u_hat_m;
          item.mode_m = objective.u_hat_m;
          Eigen::MatrixXd precision(objective.hessian_random_m);
          precision = 0.5 * (precision + precision.transpose()).eval();
          Eigen::LLT<Eigen::MatrixXd> llt;
          double jitter = 0.0;
          for (int attempt = 0; attempt < 8; ++attempt) {
            Eigen::MatrixXd regularized = precision;
            if (jitter > 0.0)
              regularized.diagonal().array() += jitter;
            llt.compute(regularized);
            if (llt.info() == Eigen::Success)
              break;
            jitter = jitter == 0.0 ? 1e-10 : jitter * 10.0;
          }
          if (llt.info() != Eigen::Success) {
            result.push_back(std::move(item));
            continue;
          }
          Eigen::VectorXd z(static_cast<Eigen::Index>(item.mode_m.size()));
          for (Eigen::Index j = 0; j < z.size(); ++j)
            z[j] = normal(rng);
          const Eigen::VectorXd conditional_deviation = llt.matrixU().solve(z);
          item.random_m = item.mode_m;
          for (size_t j = 0; j < item.random_m.size(); ++j)
            item.random_m[j] +=
                conditional_deviation[static_cast<Eigen::Index>(j)];
          item.full_m =
              quadra::merge_parameters(fixed, item.random_m, partition);
          item.valid_m = true;
        } catch (const std::exception &) {
          item.valid_m = false;
        }
        result.push_back(std::move(item));
      }
      return result;
    }));
  }
  for (size_t chain = 0; chain < futures.size(); ++chain)
    reconstructed[chain] = futures[chain].get();

  const auto csv_number = [](double value) {
    if (!std::isfinite(value))
      return std::string();
    std::ostringstream text;
    text << std::setprecision(17) << value;
    return text.str();
  };
  std::ostringstream random_draws;
  random_draws << "chain,iteration,valid";
  for (const size_t index : partition.random_indices_m)
    random_draws << ",mode_" << all_names[index] << ",draw_"
                 << all_names[index];
  random_draws << "\n";
  size_t valid_reconstructions = 0;
  for (size_t chain = 0; chain < reconstructed.size(); ++chain)
    for (size_t iteration = 0; iteration < reconstructed[chain].size();
         ++iteration) {
      const auto &item = reconstructed[chain][iteration];
      random_draws << chain + 1 << "," << iteration + 1 << ","
                   << (item.valid_m ? 1 : 0);
      for (size_t j = 0; j < partition.random_indices_m.size(); ++j)
        random_draws << "," << (item.valid_m ? csv_number(item.mode_m[j]) : "")
                     << ","
                     << (item.valid_m ? csv_number(item.random_m[j]) : "");
      random_draws << "\n";
      valid_reconstructions += item.valid_m ? 1 : 0;
    }
  quadra::write_text_file(out_dir + "/posterior_random_effect_draws.csv",
                          random_draws.str());

  const int maximum_management_draws =
      std::max(1, configured_int("QUADRA_TUNA_NUTS_MANAGEMENT_DRAWS",
                                 "sampling.management_draws", 50));
  const size_t total_draws = workflow.total_draws();
  const size_t stride = std::max<size_t>(
      1, (total_draws + static_cast<size_t>(maximum_management_draws) - 1) /
             static_cast<size_t>(maximum_management_draws));
  std::ostringstream reference_draws;
  reference_draws << "chain,iteration,valid,B0,B_MSY,MSY,F_MSY_multiplier,"
                     "B_terminal_over_B_MSY,F_status_quo_over_F_MSY\n";
  std::ostringstream projection_draws;
  projection_draws
      << "chain,iteration,scenario,projection_year,fishing_multiplier,"
         "spawning_biomass,depletion,retained_yield,discard_yield,"
         "total_yield\n";
  size_t global_draw = 0;
  int accepted_management_draws = 0;
  for (size_t chain = 0; chain < reconstructed.size(); ++chain)
    for (size_t iteration = 0; iteration < reconstructed[chain].size();
         ++iteration, ++global_draw) {
      if (global_draw % stride != 0 ||
          accepted_management_draws >= maximum_management_draws)
        continue;
      const auto &item = reconstructed[chain][iteration];
      if (!item.valid_m)
        continue;
      const auto summary = quadra::evaluate_at_parameters(
          data, controls, item.full_m, "posterior_draw");
      const auto reference = quadra::calculate_tuna_reference_points(
          data, controls, item.full_m, summary.ssb_terminal_m, 3.0, 60);
      reference_draws << chain + 1 << "," << iteration + 1 << ","
                      << (reference.valid_m ? 1 : 0) << ","
                      << csv_number(reference.b0_m) << ","
                      << csv_number(reference.b_msy_m) << ","
                      << csv_number(reference.msy_m) << ","
                      << csv_number(reference.f_msy_multiplier_m) << ","
                      << csv_number(reference.terminal_b_over_b_msy_m) << ","
                      << csv_number(reference.status_quo_f_over_f_msy_m)
                      << "\n";
      if (!reference.valid_m)
        continue;
      const double draw_sigma =
          std::exp(item.full_m[static_cast<size_t>(std::distance(
              all_names.begin(), std::find(all_names.begin(), all_names.end(),
                                           "log_sigma_recruit")))]);
      std::mt19937_64 projection_rng(options.sampler.seed + 100000u * chain +
                                     iteration);
      std::normal_distribution<double> projection_normal(0.0, 1.0);
      std::vector<double> recruitment_multipliers(10);
      for (double &value : recruitment_multipliers)
        value = std::exp(draw_sigma * projection_normal(projection_rng) -
                         0.5 * draw_sigma * draw_sigma);
      const std::vector<std::pair<std::string, double>> scenarios = {
          {"no_fishing", 0.0},
          {"half_status_quo", 0.5},
          {"status_quo", 1.0},
          {"F_MSY", reference.f_msy_multiplier_m}};
      for (const auto &scenario : scenarios)
        for (const auto &point : quadra::project_tuna_scenario(
                 data, controls, item.full_m, scenario.second, 10,
                 recruitment_multipliers))
          projection_draws << chain + 1 << "," << iteration + 1 << ","
                           << scenario.first << "," << point.projection_year_m
                           << "," << point.fishing_multiplier_m << ","
                           << point.spawning_biomass_m << ","
                           << point.depletion_m << "," << point.retained_yield_m
                           << "," << point.discard_yield_m << ","
                           << point.total_yield_m << "\n";
      ++accepted_management_draws;
    }
  quadra::write_text_file(out_dir + "/posterior_reference_points.csv",
                          reference_draws.str());
  quadra::write_text_file(out_dir + "/posterior_projection_draws.csv",
                          projection_draws.str());
  std::ostringstream reconstruction_summary;
  reconstruction_summary << "metric,value\nretained_draws," << total_draws
                         << "\nvalid_reconstructions," << valid_reconstructions
                         << "\nfailed_reconstructions,"
                         << total_draws - valid_reconstructions
                         << "\nmanagement_draws," << accepted_management_draws
                         << "\n";
  quadra::write_text_file(out_dir + "/posterior_reconstruction_summary.csv",
                          reconstruction_summary.str());
  std::filesystem::remove(out_dir + "/nuts_hessian_whitening.csv");
  if (!workflow.health.passed)
    std::cerr << "WARNING: marginal sampler health checks failed; "
                 "posterior results are not management-ready\n";
}

[[maybe_unused]] void
run_optional_nuts(const std::string &out_dir,
                  const quadra::TunaSpatialAssessmentData &data,
                  const quadra::TunaAssessmentControls &controls,
                  const quadra::TunaFitResult &fit) {
  const char *requested = std::getenv("QUADRA_TUNA_RUN_NUTS");
  if (requested == nullptr || std::string(requested) == "0")
    return;
  if (!fit.converged_m)
    throw std::runtime_error("AD-NUTS requires a converged assessment mode");

  quadra::AdvancedSpatialTunaAssessmentModel model(data, controls);
  const std::vector<std::string> centered_names = model.parameter_names();
  const std::vector<size_t> random_indices = model.random_effect_indices();
  const auto sigma_it = std::find(centered_names.begin(), centered_names.end(),
                                  "log_sigma_recruit");
  if (sigma_it == centered_names.end() || random_indices.empty())
    throw std::runtime_error("could not construct noncentered NUTS target");
  const size_t sigma_index =
      static_cast<size_t>(std::distance(centered_names.begin(), sigma_it));

  std::vector<double> noncentered_mode = fit.best_parameters_m;
  const double sigma = std::exp(noncentered_mode[sigma_index]);
  for (const size_t index : random_indices)
    noncentered_mode[index] /= sigma;

  int anchor_fleet = 1;
  if (const char *value = std::getenv("QUADRA_TUNA_ANCHOR_FLEET"))
    anchor_fleet = std::stoi(value);
  const std::vector<std::string> fixed_names = {
      "logit_steepness", "log_index_q_fleet_" + std::to_string(anchor_fleet)};
  std::vector<size_t> active_indices;
  std::vector<std::string> names;
  std::vector<double> scales;
  for (size_t i = 0; i < centered_names.size(); ++i) {
    if (std::find(fixed_names.begin(), fixed_names.end(), centered_names[i]) !=
        fixed_names.end())
      continue;
    active_indices.push_back(i);
    names.push_back(centered_names[i]);
    double scale = 0.35;
    if (std::find(random_indices.begin(), random_indices.end(), i) !=
        random_indices.end())
      scale = 0.75;
    else if (centered_names[i].find("move_logit_") == 0)
      scale = 0.4;
    else if (centered_names[i].find("sel50_raw_") == 0 ||
             centered_names[i].find("retention50_raw_") == 0)
      scale = 0.5;
    else if (centered_names[i].find("log_q_") == 0 ||
             centered_names[i].find("log_index_q_") == 0)
      scale = 0.25;
    else if (centered_names[i].find("log_sigma_") == 0)
      scale = 0.25;
    scales.push_back(scale);
  }
  const std::vector<double> initial_scales = scales;
  const size_t dimension = active_indices.size();
  std::vector<std::vector<double>> identity_whitening(
      dimension, std::vector<double>(dimension, 0.0));
  for (size_t j = 0; j < dimension; ++j)
    identity_whitening[j][j] = 1.0;
  std::vector<double> local_curvature(scales.size(),
                                      std::numeric_limits<double>::quiet_NaN());
  bool calibrate_scales = true;
  if (const char *value = std::getenv("QUADRA_TUNA_NUTS_CURVATURE_SCALE"))
    calibrate_scales = std::string(value) != "0";
  if (calibrate_scales) {
    StandardizedNoncenteredTunaPosterior calibration_target{
        &model,           sigma_index, random_indices,    active_indices,
        noncentered_mode, scales,      identity_whitening};
    std::vector<double> zero(active_indices.size(), 0.0);
    quadra::sampling::AdLogDensityEvaluator<
        StandardizedNoncenteredTunaPosterior>
        evaluator(calibration_target, zero, true);
    constexpr double h = 1e-3;
    for (size_t j = 0; j < scales.size(); ++j) {
      std::vector<double> plus = zero;
      std::vector<double> minus = zero;
      plus[j] = h;
      minus[j] = -h;
      const auto plus_eval = evaluator(plus);
      const auto minus_eval = evaluator(minus);
      if (!plus_eval.finite || !minus_eval.finite)
        continue;
      const double curvature =
          -(plus_eval.gradient[j] - minus_eval.gradient[j]) / (2.0 * h);
      local_curvature[j] = curvature;
      if (std::isfinite(curvature) && curvature > 1e-8) {
        const double adjustment =
            std::clamp(1.0 / std::sqrt(curvature), 0.01, 100.0);
        scales[j] = std::clamp(scales[j] * adjustment, 1e-6, 10.0);
      }
    }
  }
  std::ostringstream scale_audit;
  scale_audit << std::setprecision(17)
              << "parameter,initial_scale,local_curvature,calibrated_scale\n";
  for (size_t j = 0; j < scales.size(); ++j)
    scale_audit << names[j] << "," << initial_scales[j] << ","
                << local_curvature[j] << "," << scales[j] << "\n";
  quadra::write_text_file(out_dir + "/nuts_coordinate_scales.csv",
                          scale_audit.str());

  // The diagonal calibration above handles marginal scale differences.
  // Whiten the remaining local correlation using the full negative
  // Hessian of the log posterior at the fitted mode. This gives NUTS an
  // approximately isotropic coordinate system before warmup begins.
  StandardizedNoncenteredTunaPosterior hessian_target{
      &model,           sigma_index, random_indices,    active_indices,
      noncentered_mode, scales,      identity_whitening};
  std::vector<double> zero(dimension, 0.0);
  quadra::sampling::AdLogDensityEvaluator<StandardizedNoncenteredTunaPosterior>
      hessian_evaluator(hessian_target, zero, true);
  Eigen::MatrixXd curvature(dimension, dimension);
  constexpr double hessian_step = 1e-3;
  for (size_t j = 0; j < dimension; ++j) {
    std::vector<double> plus = zero;
    std::vector<double> minus = zero;
    plus[j] = hessian_step;
    minus[j] = -hessian_step;
    const auto plus_eval = hessian_evaluator(plus);
    const auto minus_eval = hessian_evaluator(minus);
    if (!plus_eval.finite || !minus_eval.finite)
      throw std::runtime_error(
          "non-finite posterior evaluation during Hessian whitening");
    for (size_t i = 0; i < dimension; ++i)
      curvature(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(j)) =
          -(plus_eval.gradient[i] - minus_eval.gradient[i]) /
          (2.0 * hessian_step);
  }
  curvature = 0.5 * (curvature + curvature.transpose()).eval();
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen(curvature);
  if (eigen.info() != Eigen::Success)
    throw std::runtime_error("posterior Hessian eigendecomposition failed");
  const double maximum_eigenvalue =
      std::max(1.0, eigen.eigenvalues().maxCoeff());
  // Weak or locally non-convex directions must not be expanded without
  // bound. A one-percent spectral floor limits the largest additional
  // whitening scale to 10x relative to the best-informed direction.
  const double eigenvalue_floor = maximum_eigenvalue * 1e-2;
  Eigen::VectorXd regularized = eigen.eigenvalues();
  size_t regularized_directions = 0;
  for (Eigen::Index j = 0; j < regularized.size(); ++j) {
    if (regularized[j] < eigenvalue_floor) {
      regularized[j] = eigenvalue_floor;
      ++regularized_directions;
    }
  }
  const Eigen::MatrixXd whitening_matrix =
      eigen.eigenvectors() * regularized.array().rsqrt().matrix().asDiagonal() *
      eigen.eigenvectors().transpose();
  std::vector<std::vector<double>> whitening(dimension,
                                             std::vector<double>(dimension));
  for (size_t i = 0; i < dimension; ++i)
    for (size_t j = 0; j < dimension; ++j)
      whitening[i][j] = whitening_matrix(static_cast<Eigen::Index>(i),
                                         static_cast<Eigen::Index>(j));
  std::ostringstream hessian_audit;
  hessian_audit << std::setprecision(17)
                << "direction,raw_eigenvalue,regularized_eigenvalue\n";
  for (Eigen::Index j = 0; j < regularized.size(); ++j)
    hessian_audit << j + 1 << "," << eigen.eigenvalues()[j] << ","
                  << regularized[j] << "\n";
  hessian_audit << "regularized_direction_count,," << regularized_directions
                << "\n";
  quadra::write_text_file(out_dir + "/nuts_hessian_whitening.csv",
                          hessian_audit.str());
  const auto env_int = [](const char *name, int fallback) {
    if (const char *value = std::getenv(name))
      return std::max(1, std::stoi(value));
    return fallback;
  };
  quadra::sampling::NutsWorkflowOptions options;
  options.chains =
      static_cast<size_t>(std::max(2, env_int("QUADRA_TUNA_NUTS_CHAINS", 4)));
  options.sampler.warmup = env_int("QUADRA_TUNA_NUTS_WARMUP", 500);
  options.sampler.samples = env_int("QUADRA_TUNA_NUTS_SAMPLES", 500);
  options.sampler.max_tree_depth =
      env_int("QUADRA_TUNA_NUTS_MAX_TREE_DEPTH", 10);
  options.sampler.target_acceptance = 0.9;
  if (const char *value = std::getenv("QUADRA_TUNA_NUTS_TARGET_ACCEPTANCE"))
    options.sampler.target_acceptance = std::stod(value);
  options.sampler.adapt_diagonal_mass = true;
  if (const char *value = std::getenv("QUADRA_TUNA_NUTS_ADAPT_MASS"))
    options.sampler.adapt_diagonal_mass = std::string(value) != "0";
  options.sampler.adapt_dense_mass = options.sampler.adapt_diagonal_mass;
  if (const char *value = std::getenv("QUADRA_TUNA_NUTS_DENSE_METRIC"))
    options.sampler.adapt_dense_mass =
        options.sampler.adapt_diagonal_mass && std::string(value) != "0";
  options.sampler.seed = 20260820;
  options.initialization_seed = 20260820;
  options.health.max_rhat = 1.05;

  std::vector<double> standardized_mode(active_indices.size(), 0.0);
  auto workflow = quadra::sampling::run_nuts_workflow(
      [&model, sigma_index, random_indices, active_indices, noncentered_mode,
       scales, whitening](size_t) {
        return StandardizedNoncenteredTunaPosterior{
            &model,           sigma_index, random_indices, active_indices,
            noncentered_mode, scales,      whitening};
      },
      standardized_mode, names, options);

  // Convert standardized non-centered draws into assessment-scale active
  // parameters before serialization and management calculations.
  for (auto &chain : workflow.fit.chains)
    for (auto &draw : chain.draws) {
      std::vector<double> full = noncentered_mode;
      const std::vector<double> standardized_draw = draw;
      for (size_t j = 0; j < active_indices.size(); ++j) {
        double transformed = 0.0;
        for (size_t k = 0; k < active_indices.size(); ++k)
          transformed += whitening[j][k] * standardized_draw[k];
        full[active_indices[j]] += scales[j] * transformed;
      }
      const double draw_sigma = std::exp(full[sigma_index]);
      for (const size_t index : random_indices)
        full[index] *= draw_sigma;
      for (size_t j = 0; j < active_indices.size(); ++j)
        draw[j] = full[active_indices[j]];
    }
  quadra::sampling::write_posterior_draws_csv(out_dir + "/posterior_draws.csv",
                                              workflow);
  quadra::sampling::write_parameter_diagnostics_csv(
      out_dir + "/posterior_parameter_diagnostics.csv", workflow);
  quadra::sampling::write_chain_diagnostics_csv(
      out_dir + "/posterior_chain_diagnostics.csv", workflow);
  quadra::sampling::write_nuts_summary_csv(out_dir + "/sampler_summary.csv",
                                           workflow);
  quadra::sampling::write_nuts_summary_csv(out_dir + "/nuts_summary.csv",
                                           workflow);

  const int maximum_management_draws =
      env_int("QUADRA_TUNA_NUTS_MANAGEMENT_DRAWS", 50);
  const size_t total_draws = workflow.total_draws();
  const size_t stride = std::max<size_t>(
      1, (total_draws + static_cast<size_t>(maximum_management_draws) - 1) /
             static_cast<size_t>(maximum_management_draws));
  std::ostringstream reference_draws;
  reference_draws << "chain,iteration,valid,B0,B_MSY,MSY,F_MSY_multiplier,"
                     "B_terminal_over_B_MSY,F_status_quo_over_F_MSY\n";
  std::ostringstream projection_draws;
  projection_draws
      << "chain,iteration,scenario,projection_year,fishing_multiplier,"
         "spawning_biomass,depletion,retained_yield,discard_yield,"
         "total_yield\n";
  const auto csv_number = [](double value) {
    if (!std::isfinite(value))
      return std::string();
    std::ostringstream text;
    text.precision(17);
    text << value;
    return text.str();
  };
  size_t global_draw = 0;
  int accepted_management_draws = 0;
  for (size_t chain = 0; chain < workflow.fit.chains.size(); ++chain) {
    const auto &draws = workflow.fit.chains[chain].draws;
    for (size_t iteration = 0; iteration < draws.size();
         ++iteration, ++global_draw) {
      if (global_draw % stride != 0 ||
          accepted_management_draws >= maximum_management_draws)
        continue;
      std::vector<double> centered = fit.best_parameters_m;
      for (size_t j = 0; j < active_indices.size(); ++j)
        centered[active_indices[j]] = draws[iteration][j];
      const double draw_sigma = std::exp(centered[sigma_index]);
      const quadra::TunaAssessmentRunSummary summary =
          quadra::evaluate_at_parameters(data, controls, centered,
                                         "posterior_draw");
      const quadra::TunaReferencePoints reference =
          quadra::calculate_tuna_reference_points(
              data, controls, centered, summary.ssb_terminal_m, 3.0, 60);
      reference_draws << chain + 1 << "," << iteration + 1 << ","
                      << (reference.valid_m ? 1 : 0) << ","
                      << csv_number(reference.b0_m) << ","
                      << csv_number(reference.b_msy_m) << ","
                      << csv_number(reference.msy_m) << ","
                      << csv_number(reference.f_msy_multiplier_m) << ","
                      << csv_number(reference.terminal_b_over_b_msy_m) << ","
                      << csv_number(reference.status_quo_f_over_f_msy_m)
                      << "\n";
      if (!reference.valid_m)
        continue;

      std::mt19937_64 projection_rng(20260820u + 100000u * chain + iteration);
      std::normal_distribution<double> normal(0.0, 1.0);
      std::vector<double> recruitment_multipliers(10);
      for (double &value : recruitment_multipliers)
        value = std::exp(draw_sigma * normal(projection_rng) -
                         0.5 * draw_sigma * draw_sigma);
      const std::vector<std::pair<std::string, double>> scenarios = {
          {"no_fishing", 0.0},
          {"half_status_quo", 0.5},
          {"status_quo", 1.0},
          {"F_MSY", reference.f_msy_multiplier_m}};
      for (const auto &scenario : scenarios) {
        const auto points = quadra::project_tuna_scenario(
            data, controls, centered, scenario.second, 10,
            recruitment_multipliers);
        for (const auto &point : points)
          projection_draws << chain + 1 << "," << iteration + 1 << ","
                           << scenario.first << "," << point.projection_year_m
                           << "," << point.fishing_multiplier_m << ","
                           << point.spawning_biomass_m << ","
                           << point.depletion_m << "," << point.retained_yield_m
                           << "," << point.discard_yield_m << ","
                           << point.total_yield_m << "\n";
      }
      ++accepted_management_draws;
    }
  }
  quadra::write_text_file(out_dir + "/posterior_reference_points.csv",
                          reference_draws.str());
  quadra::write_text_file(out_dir + "/posterior_projection_draws.csv",
                          projection_draws.str());
  if (!workflow.health.passed)
    std::cerr << "WARNING: AD-NUTS sampler health checks failed; "
                 "posterior results are not management-ready\n";
}

void write_reference_points_and_projections(
    const std::string &out_dir, const quadra::TunaSpatialAssessmentData &data,
    const quadra::TunaAssessmentControls &controls,
    const quadra::TunaFitResult &fit) {
  quadra::TunaReferencePoints reference_points;
  if (!fit.converged_m) {
    reference_points.message_m =
        "not calculated because the assessment fit did not converge";
  } else {
    reference_points = quadra::calculate_tuna_reference_points(
        data, controls, fit.best_parameters_m,
        fit.summary_m.ssb_by_year_m.back());
  }
  quadra::write_text_file(out_dir + "/reference_points.csv",
                          quadra::tuna_reference_points_csv(reference_points));

  std::ostringstream metadata;
  metadata << "valid,grid_boundary,fishing_pattern_source,message\n"
           << (reference_points.valid_m ? 1 : 0) << ","
           << (reference_points.grid_boundary_m ? 1 : 0) << ",\""
           << reference_points.fishing_pattern_source_m << "\",\""
           << reference_points.message_m << "\"\n";
  quadra::write_text_file(out_dir + "/reference_point_metadata.csv",
                          metadata.str());

  std::ostringstream projections;
  projections << "scenario,projection_year,fishing_multiplier,"
                 "spawning_biomass,depletion,retained_yield,"
                 "discard_yield,total_yield\n";
  if (reference_points.valid_m) {
    const std::vector<std::pair<std::string, double>> scenarios = {
        {"no_fishing", 0.0},
        {"half_status_quo", 0.5},
        {"status_quo", 1.0},
        {"F_MSY", reference_points.f_msy_multiplier_m}};
    for (const auto &scenario : scenarios) {
      const auto points = quadra::project_tuna_scenario(
          data, controls, fit.best_parameters_m, scenario.second, 10);
      const std::string csv =
          quadra::tuna_projection_csv(scenario.first, points);
      const size_t first_newline = csv.find('\n');
      if (first_newline != std::string::npos)
        projections << csv.substr(first_newline + 1);
    }
  }
  quadra::write_text_file(out_dir + "/projection_summary.csv",
                          projections.str());
}
} // namespace

int main(int argc, char **argv) {
  for (int i = 1; i < argc; ++i) {
    const std::string argument = argv[i];
    if (argument == "--config") {
      if (++i >= argc)
        throw std::invalid_argument("--config requires a path");
      g_config.load(argv[i]);
    } else
      throw std::invalid_argument("unknown command-line argument: " + argument);
  }
  quadra::TunaSpatialAssessmentData data;
  data.n_years_m = 12;
  data.n_ages_m = 6;
  data.n_fleets_m = 2;
  data.n_regions_m = 2;
  data.n_seasons_m = 4;
  data.spawning_fraction_m = 0.5;

  data.natural_mortality_at_age_m = {0.8, 0.5, 0.35, 0.25, 0.2, 0.18};
  data.maturity_at_age_m = {0.0, 0.2, 0.55, 0.85, 0.95, 0.98};

  data.weight_at_age_m.resize(
      static_cast<size_t>(data.n_years_m * data.n_ages_m));
  for (int y = 0; y < data.n_years_m; ++y) {
    for (int a = 0; a < data.n_ages_m; ++a) {
      data.weight_at_age_m[data.year_age_index(y, a)] =
          0.8 + 0.5 * static_cast<double>(a) + 0.02 * static_cast<double>(y);
    }
  }

  data.regional_recruit_proportions_m = {0.6, 0.4};

  data.movement_matrix_m.assign(
      static_cast<size_t>(data.n_seasons_m * data.n_regions_m *
                          data.n_regions_m),
      0.0);
  for (int s = 0; s < data.n_seasons_m; ++s) {
    data.movement_matrix_m[data.season_region_region_index(s, 0, 0)] = 0.88;
    data.movement_matrix_m[data.season_region_region_index(s, 0, 1)] = 0.12;
    data.movement_matrix_m[data.season_region_region_index(s, 1, 0)] = 0.10;
    data.movement_matrix_m[data.season_region_region_index(s, 1, 1)] = 0.90;
  }

  data.availability_surface_m.assign(static_cast<size_t>(data.n_fleets_m) *
                                         static_cast<size_t>(data.n_seasons_m) *
                                         static_cast<size_t>(data.n_regions_m) *
                                         static_cast<size_t>(data.n_ages_m),
                                     1.0);
  for (int f = 0; f < data.n_fleets_m; ++f) {
    for (int s = 0; s < data.n_seasons_m; ++s) {
      for (int r = 0; r < data.n_regions_m; ++r) {
        for (int a = 0; a < data.n_ages_m; ++a) {
          const size_t idx = data.fleet_season_region_age_index(f, s, r, a);
          const double fleet_mult = 1.0 + 0.04 * static_cast<double>(f);
          const double season_mult = 0.9 + 0.05 * static_cast<double>(s);
          const double region_mult = 0.92 + 0.08 * static_cast<double>(r);
          const double age_mult = 0.8 + 0.05 * static_cast<double>(a);
          data.availability_surface_m[idx] =
              fleet_mult * season_mult * region_mult * age_mult;
        }
      }
    }
  }

  const size_t n_fysr = static_cast<size_t>(data.n_fleets_m) *
                        static_cast<size_t>(data.n_years_m) *
                        static_cast<size_t>(data.n_seasons_m) *
                        static_cast<size_t>(data.n_regions_m);
  data.effort_m.resize(n_fysr);
  data.observed_index_m.resize(n_fysr);
  data.observed_retained_biomass_m.resize(n_fysr);
  data.observed_discard_biomass_m.resize(n_fysr);

  for (int f = 0; f < data.n_fleets_m; ++f) {
    for (int y = 0; y < data.n_years_m; ++y) {
      for (int s = 0; s < data.n_seasons_m; ++s) {
        for (int r = 0; r < data.n_regions_m; ++r) {
          const size_t idx = data.fleet_year_season_region_index(f, y, s, r);
          const double ff = static_cast<double>(f + 1);
          const double yy = static_cast<double>(y + 1);
          const double ss = static_cast<double>(s + 1);
          const double rr = static_cast<double>(r + 1);

          data.effort_m[idx] = 2.0 + 0.6 * ff + 0.2 * yy + 0.05 * ss + 0.1 * rr;
          data.observed_index_m[idx] =
              2000.0 + 100.0 * ff + 40.0 * yy + 20.0 * ss + 30.0 * rr;
          data.observed_retained_biomass_m[idx] =
              150.0 + 15.0 * ff + 8.0 * yy + 3.0 * ss + 5.0 * rr;
          data.observed_discard_biomass_m[idx] =
              20.0 + 2.5 * ff + 1.2 * yy + 0.7 * ss + 0.8 * rr;
        }
      }
    }
  }

  const size_t n_fysra = n_fysr * static_cast<size_t>(data.n_ages_m);
  data.observed_catch_numbers_m.resize(n_fysra);

  for (int f = 0; f < data.n_fleets_m; ++f) {
    for (int y = 0; y < data.n_years_m; ++y) {
      for (int s = 0; s < data.n_seasons_m; ++s) {
        for (int r = 0; r < data.n_regions_m; ++r) {
          for (int a = 0; a < data.n_ages_m; ++a) {
            const size_t idx =
                data.fleet_year_season_region_age_index(f, y, s, r, a);
            data.observed_catch_numbers_m[idx] =
                12 + 2 * a + 2 * f + y + (data.n_regions_m - r);
          }
        }
      }
    }
  }

  try {
    quadra::TunaAssessmentControls controls;
    controls.phase_m = quadra::TunaAssessmentPhase::Full;
    controls.use_priors_m =
        configured_bool("QUADRA_TUNA_USE_PRIORS", "model.use_priors", true);
    controls.movement_smoothing_weight_m = configured_double(
        "QUADRA_TUNA_MOVEMENT_SMOOTHING", "model.movement_smoothing", 2.0);
    controls.availability_smoothing_weight_m =
        configured_double("QUADRA_TUNA_AVAILABILITY_SMOOTHING",
                          "model.availability_smoothing", 1.5);
    // The supplied synthetic availability surface already contains fleet
    // scaling.  A second fleet-only multiplier is exactly confounded with
    // fishing catchability, so do not estimate both in this example.
    controls.estimate_availability_scales_m =
        configured_bool("QUADRA_TUNA_ESTIMATE_AVAILABILITY_SCALES",
                        "model.estimate_availability_scales", false);
    controls.availability_by_fleet_only_m =
        configured_bool("QUADRA_TUNA_AVAILABILITY_BY_FLEET_ONLY",
                        "model.availability_by_fleet_only", true);
    controls.share_movement_across_seasons_m =
        configured_bool("QUADRA_TUNA_SHARE_MOVEMENT",
                        "model.share_movement_across_seasons", true);
    const double observation_sd_multiplier =
        configured_double("QUADRA_TUNA_OBSERVATION_SD_MULTIPLIER",
                          "simulation.observation_sd_multiplier", 1.0);
    const double retention_prior_sd_multiplier =
        configured_double("QUADRA_TUNA_RETENTION_PRIOR_SD_MULTIPLIER",
                          "model.retention_prior_sd_multiplier", 1.0);
    if (!(observation_sd_multiplier > 0.0) ||
        !(retention_prior_sd_multiplier > 0.0)) {
      throw std::invalid_argument(
          "simulation error and retention-prior multipliers must be positive");
    }
    controls.sd_prior_retention50_raw_m *= retention_prior_sd_multiplier;
    controls.sd_prior_log_retention_slope_m *= retention_prior_sd_multiplier;
    if (const char *value = std::getenv("QUADRA_TUNA_DISABLE_COMPOSITION")) {
      controls.use_catch_composition_likelihood_m = std::string(value) == "0";
    } else
      controls.use_catch_composition_likelihood_m = configured_bool(
          "QUADRA_TUNA_USE_COMPOSITION", "model.use_composition", true);

    const bool use_model_consistent_data = configured_bool(
        "QUADRA_TUNA_MODEL_CONSISTENT_DATA", "data.model_consistent", true);
    std::vector<double> generating_parameters;
    constexpr int composition_sample_size = 200;
    const double composition_effective_n =
        configured_double("QUADRA_TUNA_COMPOSITION_EFFECTIVE_N",
                          "simulation.composition_effective_n", 8.09);
    if (!(composition_effective_n > 1.0) ||
        !(composition_effective_n < composition_sample_size)) {
      throw std::invalid_argument(
          "composition effective sample size must be between 1 and 200");
    }
    if (use_model_consistent_data) {
      quadra::TunaAssessmentControls truth_controls = controls;
      truth_controls.use_catch_composition_likelihood_m = true;
      quadra::AdvancedSpatialTunaAssessmentModel truth_model(data,
                                                             truth_controls);
      std::mt19937_64 simulation_rng(20260820);
      std::normal_distribution<double> standard_normal(0.0, 1.0);
      generating_parameters = truth_model.parameter_set().initials();
      const std::vector<std::string> truth_parameter_names =
          truth_model.parameter_set().names();
      const double recruitment_sigma = std::exp(-1.0);
      for (size_t j = 0; j < truth_parameter_names.size(); ++j) {
        if (truth_parameter_names[j].find("log_sigma_index_fleet_") == 0 ||
            truth_parameter_names[j].find("log_sigma_retained_bio_fleet_") ==
                0 ||
            truth_parameter_names[j].find("log_sigma_discard_bio_fleet_") ==
                0) {
          generating_parameters[j] += std::log(observation_sd_multiplier);
        }
        if (truth_parameter_names[j].find("recruit_dev_year_") == 0) {
          generating_parameters[j] =
              recruitment_sigma * standard_normal(simulation_rng);
        }
      }
      const double composition_alpha0 =
          composition_sample_size * (composition_effective_n - 1.0) /
          (composition_sample_size - composition_effective_n);
      const double composition_theta = 1.0 / composition_alpha0;
      for (size_t j = 0; j < truth_parameter_names.size(); ++j) {
        if (truth_parameter_names[j].find("log_theta_comp_fleet_") == 0)
          generating_parameters[j] = std::log(composition_theta);
      }
      const quadra::TunaDiagnosticBundle truth_diagnostics =
          quadra::evaluate_diagnostic_bundle(data, truth_controls,
                                             generating_parameters);
      std::vector<double> composition_probabilities(
          data.observed_catch_numbers_m.size(), 0.0);
      for (const auto &row : truth_diagnostics.observations_m) {
        const int f = row.fleet_m - 1;
        const int y = row.year_m - 1;
        const int s = row.season_m - 1;
        const int r = row.region_m - 1;
        if (row.component_m == "composition") {
          const int a = row.age_m - 1;
          const size_t idx =
              data.fleet_year_season_region_age_index(f, y, s, r, a);
          composition_probabilities[idx] = row.predicted_m;
        } else {
          const size_t idx = data.fleet_year_season_region_index(f, y, s, r);
          const double value =
              std::max(1e-8, row.predicted_m *
                                 std::exp(row.standard_deviation_m *
                                          standard_normal(simulation_rng)));
          if (row.component_m == "index")
            data.observed_index_m[idx] = value;
          else if (row.component_m == "retained")
            data.observed_retained_biomass_m[idx] = value;
          else if (row.component_m == "discard")
            data.observed_discard_biomass_m[idx] = value;
        }
      }

      for (int f = 0; f < data.n_fleets_m; ++f) {
        for (int y = 0; y < data.n_years_m; ++y) {
          for (int s = 0; s < data.n_seasons_m; ++s) {
            for (int r = 0; r < data.n_regions_m; ++r) {
              std::vector<double> probabilities(
                  static_cast<size_t>(data.n_ages_m), 0.0);
              double probability_sum = 0.0;
              const double theta = composition_theta;
              for (int a = 0; a < data.n_ages_m; ++a) {
                const size_t idx =
                    data.fleet_year_season_region_age_index(f, y, s, r, a);
                const double shape =
                    std::max(1e-8, composition_probabilities[idx] / theta);
                std::gamma_distribution<double> gamma(shape, 1.0);
                probabilities[static_cast<size_t>(a)] = gamma(simulation_rng);
                probability_sum += probabilities[static_cast<size_t>(a)];
              }
              for (double &probability : probabilities)
                probability /= probability_sum;
              std::discrete_distribution<int> draw_age(probabilities.begin(),
                                                       probabilities.end());
              for (int a = 0; a < data.n_ages_m; ++a)
                data.observed_catch_numbers_m
                    [data.fleet_year_season_region_age_index(f, y, s, r, a)] =
                    0;
              for (int draw = 0; draw < composition_sample_size; ++draw) {
                const int a = draw_age(simulation_rng);
                ++data.observed_catch_numbers_m
                      [data.fleet_year_season_region_age_index(f, y, s, r, a)];
              }
            }
          }
        }
      }
    }

    data.validate();

    quadra::TunaFitOptions fit_options;
    fit_options.multistart_m = std::max(
        1, configured_int("QUADRA_TUNA_MULTISTART", "fit.multistart", 4));
    fit_options.max_iterations_per_phase_m =
        std::max(1, configured_int("QUADRA_TUNA_MAX_PHASE_ITERATIONS",
                                   "fit.max_phase_iterations", 35));
    fit_options.initial_step_m =
        configured_double("QUADRA_TUNA_INITIAL_STEP", "fit.initial_step", 0.15);
    fit_options.min_step_m =
        configured_double("QUADRA_TUNA_MIN_STEP", "fit.min_step", 1e-8);
    fit_options.seed_m =
        configured_int("QUADRA_TUNA_FIT_SEED", "fit.seed", 20260724);
    const int anchor_fleet =
        configured_int("QUADRA_TUNA_ANCHOR_FLEET", "model.anchor_fleet", 1);
    if (anchor_fleet < 1 || anchor_fleet > data.n_fleets_m)
      throw std::invalid_argument("anchor fleet is outside the fleet range");
    fit_options.fixed_parameter_names_m = {
        "logit_steepness", "log_index_q_fleet_" + std::to_string(anchor_fleet)};
    fit_options.fixed_parameter_values_m = {0.0, -8.0};
    if (use_model_consistent_data) {
      // Match Quadra's catch-at-age convergence test: begin the public
      // exact optimizer at the supplied generating values.
      fit_options.phase_sequence_m = {quadra::TunaAssessmentPhase::Full};
    }
    fit_options.hdot_workers_m = std::max(
        0, configured_int("QUADRA_TUNA_HDOT_WORKERS", "fit.hdot_workers", 0));
    fit_options.lbfgs_print_every_m =
        std::max(0, configured_int("QUADRA_TUNA_LBFGS_PRINT_EVERY",
                                   "fit.print_every", 1));
    if (configured_bool("QUADRA_TUNA_DIRECT_FULL", "fit.direct_full", false))
      fit_options.phase_sequence_m = {quadra::TunaAssessmentPhase::Full};

    const std::string out_dir =
        configured_string("QUADRA_TUNA_OUTPUT_DIR", "output.data_dir",
                          "build/assessment_outputs/data");
    const bool resume_fit_checkpoint = configured_bool(
        "QUADRA_TUNA_LOAD_FIT_CHECKPOINT", "fit.load_checkpoint", false);
    std::filesystem::create_directories(out_dir);
    write_effective_configuration(out_dir);
    const std::string fingerprint = assessment_fingerprint(data, controls);
    quadra::TunaFitResult fit =
        resume_fit_checkpoint
            ? load_fit_checkpoint(out_dir, data, controls, fingerprint)
            : quadra::fit_spatial_assessment(data, controls, fit_options);

    const quadra::TunaAssessmentRunSummary baseline = fit.summary_m;
    const quadra::TunaDiagnosticBundle diagnostics =
        quadra::evaluate_diagnostic_bundle(data, controls,
                                           fit.best_parameters_m);

    const bool baseline_only = configured_bool("QUADRA_TUNA_BASELINE_ONLY",
                                               "run.baseline_only", false);
    if (baseline_only) {
      std::filesystem::create_directories(out_dir);
      quadra::AdvancedSpatialTunaAssessmentModel checkpoint_model(data,
                                                                  controls);
      write_fit_checkpoint(out_dir, checkpoint_model, fit, fingerprint);
      const quadra::TunaRetrospectiveResult empty_retro;
      const std::vector<quadra::TunaAssessmentRunSummary> empty_sensitivity;
      quadra::write_text_file(
          out_dir + "/acceptance_summary.json",
          quadra::fit_result_json(fit, empty_retro, empty_sensitivity));
      quadra::write_text_file(out_dir + "/residual_summary.csv",
                              quadra::residual_summary_csv(diagnostics));
      quadra::write_text_file(
          out_dir + "/likelihood_decomposition.csv",
          quadra::likelihood_decomposition_csv(diagnostics));
      quadra::write_text_file(out_dir + "/parameter_diagnostics.csv",
                              quadra::parameter_diagnostics_csv(diagnostics));
      quadra::write_text_file(out_dir + "/biomass_trajectory.csv",
                              quadra::biomass_trajectory_csv(baseline));
      write_reference_points_and_projections(out_dir, data, controls, fit);
      run_optional_marginal_nuts(out_dir, data, controls, fit);
      if (resume_fit_checkpoint) {
        std::cout << "Loaded validated fit checkpoint; optimization "
                     "and gradient audit were skipped.\n";
        return 0;
      }
      std::ostringstream sensitivity_run;
      sensitivity_run
          << "converged,iterations,gradient_norm,nll,depletion_terminal,"
             "observation_sd_multiplier,composition_effective_n,"
             "retention_prior_sd_multiplier,anchor_fleet\n"
          << (fit.converged_m ? 1 : 0) << "," << fit.total_iterations_m << ","
          << fit.gradient_norm_m << "," << baseline.nll_m << ","
          << baseline.depletion_terminal_m << "," << observation_sd_multiplier
          << "," << composition_effective_n << ","
          << retention_prior_sd_multiplier << "," << anchor_fleet << "\n";
      quadra::write_text_file(out_dir + "/sensitivity_run_summary.csv",
                              sensitivity_run.str());

      quadra::AdvancedSpatialTunaAssessmentModel audit_model(data, controls);
      const quadra::ParameterSet audit_parameters = audit_model.parameter_set();
      const quadra::ParameterPartition audit_partition =
          quadra::partition_parameters(audit_parameters);
      const std::vector<std::string> all_parameter_names =
          audit_parameters.names();
      std::ostringstream recovery;
      recovery << "parameter,generating,estimated,error,absolute_error\n";
      if (generating_parameters.size() == fit.best_parameters_m.size()) {
        for (size_t j = 0; j < generating_parameters.size(); ++j) {
          const double error =
              fit.best_parameters_m[j] - generating_parameters[j];
          recovery << all_parameter_names[j] << "," << generating_parameters[j]
                   << "," << fit.best_parameters_m[j] << "," << error << ","
                   << std::abs(error) << "\n";
        }

        const auto value_by_name = [&](const std::string &name,
                                       const std::vector<double> &values) {
          const auto it = std::find(all_parameter_names.begin(),
                                    all_parameter_names.end(), name);
          return values[static_cast<size_t>(
              std::distance(all_parameter_names.begin(), it))];
        };
        const auto add_transformed = [&](const std::string &name, double truth,
                                         double estimate) {
          const double error = estimate - truth;
          recovery << name << "," << truth << "," << estimate << "," << error
                   << "," << std::abs(error) << "\n";
        };
        for (int fleet = 1; fleet <= data.n_fleets_m; ++fleet) {
          const std::string suffix = "_fleet_" + std::to_string(fleet);
          const double truth_sel_raw =
              value_by_name("sel50_raw" + suffix, generating_parameters);
          const double fit_sel_raw =
              value_by_name("sel50_raw" + suffix, fit.best_parameters_m);
          const double truth_ret_raw =
              value_by_name("retention50_raw" + suffix, generating_parameters);
          const double fit_ret_raw =
              value_by_name("retention50_raw" + suffix, fit.best_parameters_m);
          const auto bounded_age50 = [&](double raw) {
            return 1.0 + (data.n_ages_m - 1.0) / (1.0 + std::exp(-raw));
          };
          add_transformed("sel50" + suffix, bounded_age50(truth_sel_raw),
                          bounded_age50(fit_sel_raw));
          add_transformed("sel_slope" + suffix,
                          std::exp(value_by_name("log_sel_slope" + suffix,
                                                 generating_parameters)),
                          std::exp(value_by_name("log_sel_slope" + suffix,
                                                 fit.best_parameters_m)));
          add_transformed("retention50" + suffix, bounded_age50(truth_ret_raw),
                          bounded_age50(fit_ret_raw));
          add_transformed("retention_slope" + suffix,
                          std::exp(value_by_name("log_retention_slope" + suffix,
                                                 generating_parameters)),
                          std::exp(value_by_name("log_retention_slope" + suffix,
                                                 fit.best_parameters_m)));
        }
      }
      quadra::write_text_file(out_dir + "/parameter_recovery.csv",
                              recovery.str());
      const quadra::PartitionedVector<double> audit_split =
          quadra::split_parameters(fit.best_parameters_m, audit_partition);
      const std::vector<std::string> audit_fixed_names =
          quadra::parameter_names_by_indices(audit_parameters,
                                             audit_partition.fixed_indices_m);
      const quadra::LaplaceExactLBFGSOptions audit_lbfgs_options =
          quadra::make_exact_lbfgs_options(fit_options);
      const quadra::LaplaceExactGradientResult audit_exact =
          quadra::evaluate_laplace_exact_gradient(
              audit_model, audit_split.fixed_m, audit_split.random_m,
              audit_partition, audit_lbfgs_options.gradient_m);

      std::ostringstream gradient_audit;
      gradient_audit << "parameter,value,exact_envelope_gradient,finite_"
                        "difference_laplace_gradient,omitted_logdet_gradient,"
                        "relative_error,plus_ok,minus_ok\n";

      const quadra::TunaFitResult::PhaseDiagnostic *full_diag = nullptr;
      for (const auto &phase_diag : fit.phase_diagnostics_m) {
        if (phase_diag.phase_m == quadra::TunaAssessmentPhase::Full) {
          full_diag = &phase_diag;
        }
      }

      if (full_diag != nullptr) {
        for (const std::string &name :
             full_diag->largest_gradient_parameters_m) {
          const auto name_it = std::find(audit_fixed_names.begin(),
                                         audit_fixed_names.end(), name);
          if (name_it == audit_fixed_names.end()) {
            continue;
          }
          const size_t j = static_cast<size_t>(
              std::distance(audit_fixed_names.begin(), name_it));
          const double x = audit_split.fixed_m[j];
          const double h = 1e-5 * std::max(1.0, std::abs(x)) + 1e-6;
          std::vector<double> plus_fixed = audit_split.fixed_m;
          std::vector<double> minus_fixed = audit_split.fixed_m;
          plus_fixed[j] += h;
          minus_fixed[j] -= h;

          const std::vector<double> random_start = audit_exact.u_hat_m.empty()
                                                       ? audit_split.random_m
                                                       : audit_exact.u_hat_m;
          const quadra::LaplaceObjectiveResult plus =
              quadra::evaluate_laplace_objective(
                  audit_model, plus_fixed, random_start, audit_partition,
                  audit_lbfgs_options.gradient_m.objective_m);
          const quadra::LaplaceObjectiveResult minus =
              quadra::evaluate_laplace_objective(
                  audit_model, minus_fixed, random_start, audit_partition,
                  audit_lbfgs_options.gradient_m.objective_m);
          const bool plus_ok = plus.converged_m && plus.logdet_ok_m;
          const bool minus_ok = minus.converged_m && minus.logdet_ok_m;
          const double fd =
              plus_ok && minus_ok
                  ? (plus.laplace_objective_m - minus.laplace_objective_m) /
                        (2.0 * h)
                  : std::numeric_limits<double>::quiet_NaN();
          const double exact = audit_exact.gradient_fixed_m[j];
          const double omitted = fd - exact;
          const double relative_error =
              std::abs(fd - exact) /
              std::max({1.0, std::abs(fd), std::abs(exact)});
          gradient_audit << name << "," << x << "," << exact << "," << fd << ","
                         << omitted << "," << relative_error << ","
                         << (plus_ok ? 1 : 0) << "," << (minus_ok ? 1 : 0)
                         << "\n";
        }
      }
      quadra::write_text_file(out_dir + "/laplace_gradient_audit.csv",
                              gradient_audit.str());
      std::cout << "Baseline-only convergence probe\n"
                << "converged=" << (fit.converged_m ? "true" : "false")
                << " best_start=" << fit.best_start_index_m
                << " iterations=" << fit.total_iterations_m
                << " gradient_norm=" << fit.gradient_norm_m
                << " nll=" << baseline.nll_m << "\n";
      for (const auto &phase_diag : fit.phase_diagnostics_m) {
        if (phase_diag.phase_m != quadra::TunaAssessmentPhase::Full) {
          continue;
        }
        std::cout << "Largest Full-phase gradients:\n";
        for (size_t j = 0; j < phase_diag.largest_gradient_parameters_m.size();
             ++j) {
          std::cout << "  " << phase_diag.largest_gradient_parameters_m[j]
                    << "=" << phase_diag.largest_gradient_values_m[j] << "\n";
        }
      }
      std::cout << "\nLaplace gradient audit:\n" << gradient_audit.str();
      return 0;
    }

    quadra::TunaAssessmentControls composition_off_controls = controls;
    composition_off_controls.use_catch_composition_likelihood_m = false;
    quadra::TunaFitResult composition_off_fit = quadra::fit_spatial_assessment(
        data, composition_off_controls, fit_options);
    composition_off_fit.summary_m.label_m = "composition_off";
    const quadra::TunaDiagnosticBundle composition_off_diagnostics =
        quadra::evaluate_diagnostic_bundle(
            data, composition_off_controls,
            composition_off_fit.best_parameters_m);

    bool run_intermediate_grid = true;
    if (const char *value = std::getenv("QUADRA_TUNA_GRID_INTERMEDIATE"))
      run_intermediate_grid = std::string(value) != "0";
    const std::vector<double> composition_weights =
        run_intermediate_grid ? std::vector<double>{0.0, 0.1, 0.25, 0.5, 1.0}
                              : std::vector<double>{0.0, 1.0};
    std::vector<quadra::TunaFitResult> composition_grid_fits;
    std::vector<quadra::TunaDiagnosticBundle> composition_grid_diagnostics;
    composition_grid_fits.reserve(composition_weights.size());
    composition_grid_diagnostics.reserve(composition_weights.size());

    quadra::TunaAssessmentControls zero_weight_controls = controls;
    zero_weight_controls.composition_likelihood_weight_m = 0.0;
    composition_grid_fits.push_back(composition_off_fit);
    composition_grid_diagnostics.push_back(quadra::evaluate_diagnostic_bundle(
        data, zero_weight_controls, composition_off_fit.best_parameters_m));

    const std::vector<double> intermediate_weights =
        run_intermediate_grid ? std::vector<double>{0.1, 0.25, 0.5}
                              : std::vector<double>{};
    for (double weight : intermediate_weights) {
      quadra::TunaAssessmentControls grid_controls = controls;
      grid_controls.composition_likelihood_weight_m = weight;
      quadra::TunaFitResult grid_fit =
          quadra::fit_spatial_assessment(data, grid_controls, fit_options);
      grid_fit.summary_m.label_m =
          "composition_weight_" + std::to_string(weight);
      composition_grid_diagnostics.push_back(quadra::evaluate_diagnostic_bundle(
          data, grid_controls, grid_fit.best_parameters_m));
      composition_grid_fits.push_back(std::move(grid_fit));
    }
    composition_grid_fits.push_back(fit);
    composition_grid_diagnostics.push_back(diagnostics);

    const quadra::TunaRetrospectiveResult retro =
        quadra::run_retrospective_fit(data, controls, fit_options, 2);

    std::vector<quadra::TunaSensitivityScenario> scenarios;
    scenarios.push_back({"sensitivity_M_high", 1.15, 1.0, 1.0, 1.0});
    scenarios.push_back({"sensitivity_effort_up", 1.0, 1.1, 1.0, 1.0});
    scenarios.push_back({"sensitivity_availability_low", 1.0, 1.0, 0.9, 1.0});

    const std::vector<quadra::TunaAssessmentRunSummary> sensitivities =
        quadra::run_sensitivity_fit(data, controls, fit_options, scenarios);

    int simulation_count = 50;
    if (const char *value = std::getenv("QUADRA_TUNA_SIMULATIONS")) {
      simulation_count = std::max(1, std::stoi(value));
    }
    const quadra::TunaSimulationResult sim =
        quadra::run_simulation_estimation_loop(data, controls, fit_options,
                                               simulation_count, 20260724);

    std::vector<quadra::TunaAssessmentRunSummary> management_runs;
    management_runs.push_back(baseline);
    for (const auto &run : sensitivities) {
      management_runs.push_back(run);
    }
    management_runs.push_back(composition_off_fit.summary_m);

    std::filesystem::create_directories(out_dir);
    quadra::write_text_file(out_dir + "/management_summary.csv",
                            quadra::management_summary_csv(management_runs));
    quadra::write_text_file(out_dir + "/retrospective_summary.csv",
                            quadra::retrospective_summary_csv(retro));
    quadra::write_text_file(out_dir + "/simulation_summary.csv",
                            quadra::simulation_summary_csv(sim));
    quadra::write_text_file(out_dir + "/acceptance_summary.json",
                            quadra::fit_result_json(fit, retro, sensitivities));
    quadra::write_text_file(out_dir + "/observation_diagnostics.csv",
                            quadra::observation_diagnostics_csv(diagnostics));
    quadra::write_text_file(out_dir + "/residual_summary.csv",
                            quadra::residual_summary_csv(diagnostics));
    quadra::write_text_file(
        out_dir + "/stratified_residual_summary.csv",
        quadra::stratified_residual_summary_csv(diagnostics));
    quadra::write_text_file(out_dir + "/likelihood_decomposition.csv",
                            quadra::likelihood_decomposition_csv(diagnostics));
    quadra::write_text_file(out_dir + "/parameter_diagnostics.csv",
                            quadra::parameter_diagnostics_csv(diagnostics));
    quadra::write_text_file(out_dir + "/biomass_trajectory.csv",
                            quadra::biomass_trajectory_csv(baseline));
    write_reference_points_and_projections(out_dir, data, controls, fit);
    quadra::write_text_file(out_dir + "/catchability_availability.csv",
                            quadra::catchability_availability_csv(diagnostics));
    quadra::write_text_file(
        out_dir + "/composition_off_residual_summary.csv",
        quadra::residual_summary_csv(composition_off_diagnostics));
    quadra::write_text_file(
        out_dir + "/composition_off_likelihood_decomposition.csv",
        quadra::likelihood_decomposition_csv(composition_off_diagnostics));
    quadra::write_text_file(out_dir + "/composition_weight_grid.csv",
                            quadra::composition_weight_grid_csv(
                                composition_weights, composition_grid_fits,
                                composition_grid_diagnostics, data.n_fleets_m));

    quadra::AdvancedSpatialTunaAssessmentModel model(data);
    const quadra::ParameterSet p = model.parameter_set();
    const std::vector<double> initial = p.initials();

    quadra::ModelReportContext ctx;
    const double nll = model.evaluate(initial, ctx);

    std::cout << "Initial objective (NLL): " << nll << "\n";
    for (const auto &rv : ctx.reports().values()) {
      std::cout << rv.name_m << " = " << rv.value_m;
      if (rv.requires_se_m) {
        std::cout << " (adreport)";
      }
      std::cout << "\n";
    }

    std::cout << "\nAcceptance summary\n";
    std::cout << "Baseline formulation: intermediate (fleet-only availability, "
                 "shared movement)\n";
    std::cout << "Best-start index: " << fit.best_start_index_m << "\n";
    std::cout << "Objective evaluations: " << fit.objective_evaluations_m
              << "\n";
    std::cout << "Baseline depletion: " << baseline.depletion_terminal_m
              << "\n";
    std::cout << "Retrospective Mohn's rho (SSB): " << retro.mohns_rho_ssb_m
              << "\n";

    for (const auto &point : retro.points_m) {
      std::cout << "  peel=" << point.peel_m
                << " years=" << point.summary_m.n_years_m
                << " depletion=" << point.summary_m.depletion_terminal_m
                << " nll=" << point.summary_m.nll_m << "\n";
    }

    std::cout << "\nManagement summary CSV\n";
    std::cout << quadra::management_summary_csv(management_runs);
    std::cout << "\nSimulation mean depletion bias: "
              << sim.mean_depletion_bias_m << "\n";
    std::cout << "Simulation median depletion bias: "
              << sim.median_depletion_bias_m << "\n";
    std::cout << "Simulation depletion bias p10/p90: "
              << sim.p10_depletion_bias_m << " / " << sim.p90_depletion_bias_m
              << "\n";
    std::cout << "Low depletion failure rate (< "
              << sim.low_depletion_threshold_m
              << "): " << sim.low_depletion_rate_m << " ("
              << sim.low_depletion_count_m << "/" << sim.n_finite_bias_m
              << ")\n";

    std::cout << "\nResidual diagnostics by component and fleet\n";
    std::cout << quadra::residual_summary_csv(diagnostics);
    std::cout << "\nLikelihood decomposition\n";
    std::cout << quadra::likelihood_decomposition_csv(diagnostics);
    std::cout << "\nComposition weight grid\n";
    std::cout << quadra::composition_weight_grid_csv(
        composition_weights, composition_grid_fits,
        composition_grid_diagnostics, data.n_fleets_m);

    std::cout << "\nWrote outputs to: " << out_dir << "\n";
  } catch (const std::exception &ex) {
    std::cerr << "Error: " << ex.what() << "\n";
    return 1;
  }

  return 0;
}
