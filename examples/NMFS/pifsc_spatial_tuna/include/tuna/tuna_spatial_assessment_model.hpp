#ifndef QUADRA_TUNA_SPATIAL_ASSESSMENT_MODEL_HPP
#define QUADRA_TUNA_SPATIAL_ASSESSMENT_MODEL_HPP
#pragma once

#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "core/model/parameter.hpp"
#include "core/model/quadra_model.hpp"
#include "math/distributions.hpp"
#include "tuna_spatial_data.hpp"

// Compatibility for Quadra's third-order exact-Laplace path.  The current
// distribution layer supplies lgamma for AReal but not ThirdOrderScalar.
namespace had {
inline Real tetragamma_approx(Real x) {
  Real result = Real(0.0);
  while (x < Real(6.0)) {
    result -= Real(2.0) / (x * x * x);
    x += Real(1.0);
  }
  const Real inv = Real(1.0) / x;
  const Real inv2 = inv * inv;
  const Real inv3 = inv2 * inv;
  const Real inv4 = inv2 * inv2;
  const Real inv6 = inv3 * inv3;
  const Real inv8 = inv4 * inv4;
  result += -inv2 - inv3 - Real(0.5) * inv4 + (Real(1.0) / Real(6.0)) * inv6 +
            (Real(1.0) / Real(6.0)) * inv8;
  return result;
}

inline ThirdOrderScalar lgamma(const ThirdOrderScalar &x) {
  return unary_chain(x, std::lgamma(x.val), ::digamma(x.val), ::trigamma(x.val),
                     tetragamma_approx(x.val));
}
} // namespace had

namespace quadra {

enum class TunaAssessmentPhase {
  InitializeRecruitment,
  InitializeCatchability,
  InitializeMovement,
  Full
};

struct TunaAssessmentControls {
  TunaAssessmentPhase phase_m = TunaAssessmentPhase::Full;
  bool estimate_availability_scales_m = true;
  bool availability_by_fleet_only_m = false;
  bool share_movement_across_seasons_m = false;
  bool use_index_likelihood_m = true;
  bool use_catch_composition_likelihood_m = true;
  bool use_retained_biomass_likelihood_m = true;
  bool use_discard_biomass_likelihood_m = true;
  bool use_catch_conditioning_m = false;
  bool use_priors_m = true;
  bool report_observation_predictions_m = false;
  double composition_likelihood_weight_m = 1.0;

  double sd_prior_log_q_m = 2.0;
  double sd_prior_log_sel50_m = 1.5;
  double sd_prior_log_sel_slope_m = 1.2;
  double sd_prior_retention50_raw_m = 1.0;
  double sd_prior_log_retention_slope_m = 0.7;
  double sd_prior_log_availability_scale_m = 1.5;
  double sd_prior_move_logit_m = 2.0;
  double sd_prior_log_sigma_m = 1.5;

  double movement_smoothing_weight_m = 1.0;
  double availability_smoothing_weight_m = 1.0;
};

class AdvancedSpatialTunaAssessmentModel
    : public QuadraModel<AdvancedSpatialTunaAssessmentModel> {
public:
  explicit AdvancedSpatialTunaAssessmentModel(
      TunaSpatialAssessmentData data,
      TunaAssessmentControls controls = TunaAssessmentControls())
      : data_m(std::move(data)), controls_m(controls) {
    data_m.validate();
    if (!std::isfinite(controls_m.composition_likelihood_weight_m) ||
        controls_m.composition_likelihood_weight_m < 0.0) {
      throw std::invalid_argument(
          "AdvancedSpatialTunaAssessmentModel: composition likelihood weight "
          "must be finite and non-negative");
    }
  }

  void set_controls(const TunaAssessmentControls &controls) {
    controls_m = controls;
  }

  const TunaAssessmentControls &controls() const { return controls_m; }

  ParameterSet parameter_set() const {
    ParameterSet p;
    p.add("log_r0", 10.0);
    p.add("logit_steepness", 0.0);

    for (int f = 0; f < data_m.n_fleets_m; ++f) {
      p.add("log_q_fleet_" + std::to_string(f + 1), -8.0);
      p.add("log_index_q_fleet_" + std::to_string(f + 1), -8.0);
      p.add("sel50_raw_fleet_" + std::to_string(f + 1), 0.0);
      p.add("log_sel_slope_fleet_" + std::to_string(f + 1), 0.0);
      p.add("retention50_raw_fleet_" + std::to_string(f + 1), 0.0);
      p.add("log_retention_slope_fleet_" + std::to_string(f + 1), -0.2);
      p.add("log_sigma_index_fleet_" + std::to_string(f + 1), -1.0);
      p.add("log_theta_comp_fleet_" + std::to_string(f + 1), -2.0);
      p.add("log_sigma_retained_bio_fleet_" + std::to_string(f + 1), -0.8);
      p.add("log_sigma_discard_bio_fleet_" + std::to_string(f + 1), -0.8);

      if (controls_m.estimate_availability_scales_m) {
        if (controls_m.availability_by_fleet_only_m) {
          p.add("log_availability_scale_fleet_" + std::to_string(f + 1), 0.0);
        } else {
          for (int s = 0; s < data_m.n_seasons_m; ++s) {
            for (int r = 0; r < data_m.n_regions_m; ++r) {
              p.add("log_availability_scale_fleet_" + std::to_string(f + 1) +
                        "_season_" + std::to_string(s + 1) + "_region_" +
                        std::to_string(r + 1),
                    0.0);
            }
          }
        }
      }
    }

    const int movement_seasons = movement_parameter_season_count();
    for (int s = 0; s < movement_seasons; ++s) {
      for (int from = 0; from < data_m.n_regions_m; ++from) {
        const int n_free = data_m.n_regions_m - 1;
        for (int to = 0; to < n_free; ++to) {
          const double init = controls_m.share_movement_across_seasons_m
                                  ? movement_initial_logit_shared(from, to)
                                  : movement_initial_logit(s, from, to);
          p.add("move_logit_season_" + std::to_string(s + 1) + "_from_" +
                    std::to_string(from + 1) + "_to_" + std::to_string(to + 1),
                init);
        }
      }
    }

    p.add("log_sigma_recruit", -1.0);

    for (int y = 0; y < data_m.n_years_m; ++y) {
      p.add("recruit_dev_year_" + std::to_string(y + 1), 0.0,
            ParameterTransform::Identity, true);
    }

    return p;
  }

  std::vector<std::string> parameter_names_impl() const {
    return parameter_set().names();
  }

  std::vector<size_t> random_effect_indices() const {
    const size_t n_fixed = static_cast<size_t>(
        2 + per_fleet_parameter_count() * data_m.n_fleets_m +
        movement_logit_count() + 1);
    std::vector<size_t> indices;
    indices.reserve(static_cast<size_t>(data_m.n_years_m));

    for (int y = 0; y < data_m.n_years_m; ++y) {
      indices.push_back(n_fixed + static_cast<size_t>(y));
    }

    return indices;
  }

  template <typename Type>
  Type evaluate_impl(const std::vector<Type> &parameters,
                     ModelReportContext &ctx) const {
    // Read this objective in the same order as an assessment report:
    // establish what data streams are active, validate the fitted state,
    // then tell the complete population-and-observation story.
    const ObservationProcesses observations =
        observation_processes_for_current_phase();
    validate_objective_inputs(parameters);
    return evaluate_assessment_story(parameters, observations, ctx);
  }

private:
  struct ObservationProcesses {
    bool index_m;
    bool composition_m;
    bool retained_biomass_m;
    bool discard_biomass_m;
  };

  ObservationProcesses observation_processes_for_current_phase() const {
    ObservationProcesses out{controls_m.use_index_likelihood_m,
                             controls_m.use_catch_composition_likelihood_m,
                             controls_m.use_retained_biomass_likelihood_m,
                             controls_m.use_discard_biomass_likelihood_m};

    if (controls_m.phase_m == TunaAssessmentPhase::InitializeRecruitment) {
      out.composition_m = false;
      out.retained_biomass_m = false;
      out.discard_biomass_m = false;
    } else if (controls_m.phase_m ==
               TunaAssessmentPhase::InitializeCatchability) {
      out.composition_m = false;
    }
    return out;
  }

  template <typename Type>
  void validate_objective_inputs(const std::vector<Type> &parameters) const {
    if (controls_m.use_catch_conditioning_m &&
        data_m.observed_total_catch_m.empty()) {
      throw std::invalid_argument(
          "catch conditioning requested without observed total catch data");
    }

    if (parameters.size() != parameter_set().size()) {
      throw std::invalid_argument(
          "AdvancedSpatialTunaAssessmentModel: parameter vector length "
          "mismatch");
    }
  }

  template <typename Type>
  Type evaluate_assessment_story(const std::vector<Type> &parameters,
                                 const ObservationProcesses &observations,
                                 ModelReportContext &ctx) const {
    const int y_count = data_m.n_years_m;
    const int a_count = data_m.n_ages_m;
    const int f_count = data_m.n_fleets_m;
    const int r_count = data_m.n_regions_m;
    const int s_count = data_m.n_seasons_m;
    const bool use_index = observations.index_m;
    const bool use_comp = observations.composition_m;
    const bool use_retained = observations.retained_biomass_m;
    const bool use_discard = observations.discard_biomass_m;

    size_t pos = 0;
    const Type log_r0 = parameters[pos++];
    const Type logit_h = parameters[pos++];

    std::vector<Type> log_q(static_cast<size_t>(f_count));
    std::vector<Type> log_index_q(static_cast<size_t>(f_count));
    std::vector<Type> log_sel50(static_cast<size_t>(f_count));
    std::vector<Type> log_sel_slope(static_cast<size_t>(f_count));
    std::vector<Type> log_ret50(static_cast<size_t>(f_count));
    std::vector<Type> log_ret_slope(static_cast<size_t>(f_count));
    std::vector<Type> log_sigma_index(static_cast<size_t>(f_count));
    std::vector<Type> log_theta_comp(static_cast<size_t>(f_count));
    std::vector<Type> log_sigma_ret_bio(static_cast<size_t>(f_count));
    std::vector<Type> log_sigma_dis_bio(static_cast<size_t>(f_count));
    std::vector<Type> log_availability_scale(
        static_cast<size_t>(f_count * s_count * r_count));
    std::vector<Type> move_logits(
        static_cast<size_t>(s_count * r_count * (r_count - 1)));

    for (int f = 0; f < f_count; ++f) {
      const size_t ff = static_cast<size_t>(f);
      log_q[ff] = parameters[pos++];
      log_index_q[ff] = parameters[pos++];
      log_sel50[ff] = parameters[pos++];
      log_sel_slope[ff] = parameters[pos++];
      log_ret50[ff] = parameters[pos++];
      log_ret_slope[ff] = parameters[pos++];
      log_sigma_index[ff] = parameters[pos++];
      log_theta_comp[ff] = parameters[pos++];
      log_sigma_ret_bio[ff] = parameters[pos++];
      log_sigma_dis_bio[ff] = parameters[pos++];

      if (controls_m.estimate_availability_scales_m &&
          controls_m.availability_by_fleet_only_m) {
        const Type fleet_log_avail = parameters[pos++];
        for (int s = 0; s < s_count; ++s) {
          for (int r = 0; r < r_count; ++r) {
            log_availability_scale[fleet_season_region_index(f, s, r)] =
                fleet_log_avail;
          }
        }
      } else {
        for (int s = 0; s < s_count; ++s) {
          for (int r = 0; r < r_count; ++r) {
            if (controls_m.estimate_availability_scales_m) {
              log_availability_scale[fleet_season_region_index(f, s, r)] =
                  parameters[pos++];
            } else {
              log_availability_scale[fleet_season_region_index(f, s, r)] =
                  Type(0.0);
            }
          }
        }
      }
    }

    const int movement_seasons = movement_parameter_season_count();
    for (int s = 0; s < movement_seasons; ++s) {
      for (int from = 0; from < r_count; ++from) {
        for (int to = 0; to < r_count - 1; ++to) {
          move_logits[movement_logit_index(s, from, to)] = parameters[pos++];
        }
      }
    }

    const Type log_sigma_recruit = parameters[pos++];

    std::vector<Type> recruit_dev(static_cast<size_t>(y_count));
    for (int y = 0; y < y_count; ++y) {
      recruit_dev[static_cast<size_t>(y)] = parameters[pos++];
    }

    Type nll = Type(0.0);
    Type nll_prior = Type(0.0);
    Type nll_availability_penalty = Type(0.0);
    Type nll_movement_penalty = Type(0.0);
    Type nll_index = Type(0.0);
    Type nll_composition = Type(0.0);
    Type nll_retained = Type(0.0);
    Type nll_discard = Type(0.0);
    Type nll_recruitment = Type(0.0);

    const Type r0 = exp(log_r0);
    const Type steepness = Type(0.2) + Type(0.8) * logistic(logit_h);
    const Type sigma_recruit = exp(log_sigma_recruit);

    std::vector<Type> q(static_cast<size_t>(f_count));
    std::vector<Type> index_q(static_cast<size_t>(f_count));
    std::vector<Type> sel50(static_cast<size_t>(f_count));
    std::vector<Type> sel_slope(static_cast<size_t>(f_count));
    std::vector<Type> ret50(static_cast<size_t>(f_count));
    std::vector<Type> ret_slope(static_cast<size_t>(f_count));
    std::vector<Type> sigma_index(static_cast<size_t>(f_count));
    std::vector<Type> theta_comp(static_cast<size_t>(f_count));
    std::vector<Type> sigma_ret_bio(static_cast<size_t>(f_count));
    std::vector<Type> sigma_dis_bio(static_cast<size_t>(f_count));
    std::vector<Type> availability_scale(
        static_cast<size_t>(f_count * s_count * r_count));
    std::vector<Type> move_probs(
        static_cast<size_t>(s_count * r_count * r_count));

    for (int f = 0; f < f_count; ++f) {
      const size_t ff = static_cast<size_t>(f);
      q[ff] = exp(log_q[ff]);
      index_q[ff] = exp(log_index_q[ff]);
      sel50[ff] = bounded_age50(log_sel50[ff]);
      sel_slope[ff] = exp(log_sel_slope[ff]);
      ret50[ff] = bounded_age50(log_ret50[ff]);
      ret_slope[ff] = exp(log_ret_slope[ff]);
      sigma_index[ff] = exp(log_sigma_index[ff]);
      theta_comp[ff] = exp(log_theta_comp[ff]);
      sigma_ret_bio[ff] = exp(log_sigma_ret_bio[ff]);
      sigma_dis_bio[ff] = exp(log_sigma_dis_bio[ff]);

      for (int s = 0; s < s_count; ++s) {
        for (int r = 0; r < r_count; ++r) {
          availability_scale[fleet_season_region_index(f, s, r)] =
              exp(log_availability_scale[fleet_season_region_index(f, s, r)]);
        }
      }
    }

    for (int s = 0; s < s_count; ++s) {
      for (int from = 0; from < r_count; ++from) {
        std::vector<Type> numerators(static_cast<size_t>(r_count), Type(0.0));
        Type denom = Type(0.0);

        for (int to = 0; to < r_count - 1; ++to) {
          const int s_param = movement_parameter_season_index(s);
          const Type num =
              exp(move_logits[movement_logit_index(s_param, from, to)]);
          numerators[static_cast<size_t>(to)] = num;
          denom += num;
        }

        numerators[static_cast<size_t>(r_count - 1)] = Type(1.0);
        denom += Type(1.0);

        for (int to = 0; to < r_count; ++to) {
          move_probs[movement_probability_index(s, from, to)] =
              numerators[static_cast<size_t>(to)] / denom;
        }
      }
    }

    if (controls_m.use_priors_m) {
      const Type sd_log_q = Type(controls_m.sd_prior_log_q_m);
      const Type sd_log_sel50 = Type(controls_m.sd_prior_log_sel50_m);
      const Type sd_log_sel_slope = Type(controls_m.sd_prior_log_sel_slope_m);
      const Type sd_retention50 = Type(controls_m.sd_prior_retention50_raw_m);
      const Type sd_log_retention_slope =
          Type(controls_m.sd_prior_log_retention_slope_m);
      const Type sd_log_avail =
          Type(controls_m.sd_prior_log_availability_scale_m);
      const Type sd_move = Type(controls_m.sd_prior_move_logit_m);
      const Type sd_log_sigma = Type(controls_m.sd_prior_log_sigma_m);

      for (int f = 0; f < f_count; ++f) {
        const size_t ff = static_cast<size_t>(f);
        nll_prior -= dnorm(log_q[ff], Type(0.0), sd_log_q, true);
        nll -= dnorm(log_q[ff], Type(0.0), sd_log_q, true);
        nll_prior -= dnorm(log_index_q[ff], Type(0.0), sd_log_q, true);
        nll -= dnorm(log_index_q[ff], Type(0.0), sd_log_q, true);
        nll_prior -= dnorm(log_sel50[ff], Type(0.0), sd_log_sel50, true);
        nll -= dnorm(log_sel50[ff], Type(0.0), sd_log_sel50, true);
        nll_prior -=
            dnorm(log_sel_slope[ff], Type(0.0), sd_log_sel_slope, true);
        nll -= dnorm(log_sel_slope[ff], Type(0.0), sd_log_sel_slope, true);
        nll_prior -= dnorm(log_ret50[ff], Type(0.0), sd_retention50, true);
        nll -= dnorm(log_ret50[ff], Type(0.0), sd_retention50, true);
        nll_prior -=
            dnorm(log_ret_slope[ff], Type(-0.2), sd_log_retention_slope, true);
        nll -=
            dnorm(log_ret_slope[ff], Type(-0.2), sd_log_retention_slope, true);
        nll_prior -= dnorm(log_sigma_index[ff], Type(-1.0), sd_log_sigma, true);
        nll -= dnorm(log_sigma_index[ff], Type(-1.0), sd_log_sigma, true);
        nll_prior -=
            dnorm(log_sigma_ret_bio[ff], Type(-1.0), sd_log_sigma, true);
        nll -= dnorm(log_sigma_ret_bio[ff], Type(-1.0), sd_log_sigma, true);
        nll_prior -=
            dnorm(log_sigma_dis_bio[ff], Type(-1.0), sd_log_sigma, true);
        nll -= dnorm(log_sigma_dis_bio[ff], Type(-1.0), sd_log_sigma, true);

        if (controls_m.estimate_availability_scales_m &&
            controls_m.availability_by_fleet_only_m) {
          const Type x =
              log_availability_scale[fleet_season_region_index(f, 0, 0)];
          nll_prior -= dnorm(x, Type(0.0), sd_log_avail, true);
          nll -= dnorm(x, Type(0.0), sd_log_avail, true);
        } else {
          for (int s = 0; s < s_count; ++s) {
            for (int r = 0; r < r_count; ++r) {
              const Type x =
                  log_availability_scale[fleet_season_region_index(f, s, r)];
              nll_prior -= dnorm(x, Type(0.0), sd_log_avail, true);
              nll -= dnorm(x, Type(0.0), sd_log_avail, true);
            }
          }
        }
      }

      for (int s = 0; s < s_count; ++s) {
        for (int from = 0; from < r_count; ++from) {
          for (int to = 0; to < r_count - 1; ++to) {
            const int s_param = movement_parameter_season_index(s);
            const Type x = move_logits[movement_logit_index(s_param, from, to)];
            nll_prior -= dnorm(x, Type(0.0), sd_move, true);
            nll -= dnorm(x, Type(0.0), sd_move, true);
          }
        }
      }
    }

    if (controls_m.availability_smoothing_weight_m > 0.0) {
      const Type w = Type(controls_m.availability_smoothing_weight_m);
      for (int f = 0; f < f_count; ++f) {
        for (int s = 1; s < s_count; ++s) {
          for (int r = 0; r < r_count; ++r) {
            const Type x1 =
                log_availability_scale[fleet_season_region_index(f, s, r)];
            const Type x0 =
                log_availability_scale[fleet_season_region_index(f, s - 1, r)];
            const Type d = x1 - x0;
            nll_availability_penalty += Type(0.5) * w * d * d;
            nll += Type(0.5) * w * d * d;
          }
        }
      }
    }

    if (controls_m.movement_smoothing_weight_m > 0.0) {
      const Type w = Type(controls_m.movement_smoothing_weight_m);
      const int movement_smoothing_seasons = movement_parameter_season_count();
      for (int s = 1; s < movement_smoothing_seasons; ++s) {
        for (int from = 0; from < r_count; ++from) {
          for (int to = 0; to < r_count - 1; ++to) {
            const Type x1 = move_logits[movement_logit_index(s, from, to)];
            const Type x0 = move_logits[movement_logit_index(s - 1, from, to)];
            const Type d = x1 - x0;
            nll_movement_penalty += Type(0.5) * w * d * d;
            nll += Type(0.5) * w * d * d;
          }
        }
      }
    }

    std::vector<Type> state = initialize_state(r0, sigma_recruit, recruit_dev);
    const Type ssb0 = spawning_biomass(state, 0);
    std::vector<Type> ssb_by_year(static_cast<size_t>(y_count), Type(0.0));

    std::vector<Type> selectivity(static_cast<size_t>(f_count) *
                                  static_cast<size_t>(a_count));
    std::vector<Type> retention(static_cast<size_t>(f_count) *
                                static_cast<size_t>(a_count));

    for (int f = 0; f < f_count; ++f) {
      for (int a = 0; a < a_count; ++a) {
        const size_t fa_idx = fleet_age_index(f, a);
        const Type age = Type(static_cast<double>(a) + 1.0);
        selectivity[fa_idx] =
            normalized_logistic_age(age, sel50[static_cast<size_t>(f)],
                                    sel_slope[static_cast<size_t>(f)]);
        retention[fa_idx] = logistic_age(age, ret50[static_cast<size_t>(f)],
                                         ret_slope[static_cast<size_t>(f)]);
      }
    }

    for (int y = 0; y < y_count; ++y) {
      ssb_by_year[static_cast<size_t>(y)] = spawning_biomass(state, y);

      for (int s = 0; s < s_count; ++s) {
        std::vector<Type> survivors(state.size(), Type(0.0));
        std::vector<Type> moved(state.size(), Type(0.0));

        std::vector<Type> pred_capture_age(static_cast<size_t>(f_count) *
                                               static_cast<size_t>(r_count) *
                                               static_cast<size_t>(a_count),
                                           Type(0.0));
        std::vector<Type> pred_vulnerable_bio(static_cast<size_t>(f_count) *
                                                  static_cast<size_t>(r_count),
                                              Type(0.0));
        std::vector<Type> pred_retained_bio(static_cast<size_t>(f_count) *
                                                static_cast<size_t>(r_count),
                                            Type(0.0));
        std::vector<Type> pred_discard_bio(static_cast<size_t>(f_count) *
                                               static_cast<size_t>(r_count),
                                           Type(0.0));
        std::vector<Type> conditioned_harvest_scale(
            static_cast<size_t>(f_count) * static_cast<size_t>(r_count),
            Type(0.0));

        if (controls_m.use_catch_conditioning_m) {
          for (int f = 0; f < f_count; ++f) {
            for (int r = 0; r < r_count; ++r) {
              Type vulnerable = Type(1e-6);
              for (int a = 0; a < a_count; ++a) {
                const size_t ra_idx = data_m.region_age_index(r, a);
                const size_t fa_idx = fleet_age_index(f, a);
                const Type avail = effective_availability(
                    f, s, r, a,
                    availability_scale[fleet_season_region_index(f, s, r)]);
                Type contribution = state[ra_idx] * selectivity[fa_idx] * avail;
                if (data_m.catch_units_m[static_cast<size_t>(f)] == 1) {
                  contribution *=
                      Type(data_m.weight_at_age_m[data_m.year_age_index(y, a)]);
                }
                vulnerable += contribution;
              }
              const size_t fysr_idx =
                  data_m.fleet_year_season_region_index(f, y, s, r);
              conditioned_harvest_scale[fleet_region_index(f, r)] =
                  Type(data_m.observed_total_catch_m[fysr_idx]) / vulnerable;
              if (controls_m.report_observation_predictions_m) {
                ctx.report("diag_conditioned_scale_y" + std::to_string(y + 1) +
                               "_s" + std::to_string(s + 1) + "_f" +
                               std::to_string(f + 1) + "_r" +
                               std::to_string(r + 1),
                           conditioned_harvest_scale[fleet_region_index(f, r)]);
              }
            }
          }
        }

        for (int r = 0; r < r_count; ++r) {
          for (int a = 0; a < a_count; ++a) {
            const size_t ra_idx = data_m.region_age_index(r, a);
            const Type n = state[ra_idx];
            const Type m =
                Type(
                    data_m.natural_mortality_at_age_m[static_cast<size_t>(a)]) /
                Type(static_cast<double>(s_count));

            Type f_total = Type(0.0);
            for (int f = 0; f < f_count; ++f) {
              const size_t fysr_idx =
                  data_m.fleet_year_season_region_index(f, y, s, r);
              const size_t fa_idx = fleet_age_index(f, a);
              const Type avail = effective_availability(
                  f, s, r, a,
                  availability_scale[fleet_season_region_index(f, s, r)]);
              const Type fleet_rate =
                  controls_m.use_catch_conditioning_m
                      ? conditioned_harvest_scale[fleet_region_index(f, r)] *
                            selectivity[fa_idx] * avail
                      : q[static_cast<size_t>(f)] *
                            Type(data_m.effort_m[fysr_idx]) *
                            selectivity[fa_idx] * avail;
              f_total += fleet_rate;
            }

            const Type z = m + f_total + Type(1e-12);
            const Type exploitation = (Type(1.0) - exp(-z)) / z;
            const Type weight =
                Type(data_m.weight_at_age_m[data_m.year_age_index(y, a)]);

            for (int f = 0; f < f_count; ++f) {
              const size_t fysr_idx =
                  data_m.fleet_year_season_region_index(f, y, s, r);
              const size_t fa_idx = fleet_age_index(f, a);
              const size_t fra_idx = fleet_region_age_index(f, r, a);
              const size_t fr_idx = fleet_region_index(f, r);
              const Type avail = effective_availability(
                  f, s, r, a,
                  availability_scale[fleet_season_region_index(f, s, r)]);

              const Type f_fleet = controls_m.use_catch_conditioning_m
                                       ? conditioned_harvest_scale[fr_idx] *
                                             selectivity[fa_idx] * avail
                                       : q[static_cast<size_t>(f)] *
                                             Type(data_m.effort_m[fysr_idx]) *
                                             selectivity[fa_idx] * avail;
              const Type capture_n = controls_m.use_catch_conditioning_m
                                         ? n * f_fleet
                                         : n * f_fleet * exploitation;
              const Type ret = retention[fa_idx];

              pred_capture_age[fra_idx] += capture_n;
              pred_vulnerable_bio[fr_idx] += n * weight * selectivity[fa_idx];
              pred_retained_bio[fr_idx] += capture_n * ret * weight;
              pred_discard_bio[fr_idx] +=
                  capture_n * (Type(1.0) - ret) * weight;
            }

            if (controls_m.use_catch_conditioning_m) {
              const Type remaining = Type(1.0) - f_total;
              const Type positive_remaining =
                  log(Type(1.0) + exp(Type(40.0) * remaining)) / Type(40.0);
              survivors[ra_idx] = n * positive_remaining * exp(-m);
            } else {
              survivors[ra_idx] = n * exp(-z);
            }
          }
        }

        for (int f = 0; f < f_count; ++f) {
          for (int r = 0; r < r_count; ++r) {
            const size_t fr_idx = fleet_region_index(f, r);
            const size_t fysr_idx =
                data_m.fleet_year_season_region_index(f, y, s, r);

            std::vector<Type> pred_comp(static_cast<size_t>(a_count),
                                        Type(0.0));
            std::vector<int> obs_comp(static_cast<size_t>(a_count), 0);

            Type pred_total = Type(0.0);
            int obs_total = 0;
            for (int a = 0; a < a_count; ++a) {
              const size_t fra_idx = fleet_region_age_index(f, r, a);
              pred_total += pred_capture_age[fra_idx];

              const int obs = data_m.observed_catch_numbers_m
                                  [data_m.fleet_year_season_region_age_index(
                                      f, y, s, r, a)];
              obs_comp[static_cast<size_t>(a)] = obs;
              obs_total += obs;
            }

            if (use_comp && obs_total > 0) {
              for (int a = 0; a < a_count; ++a) {
                const size_t fra_idx = fleet_region_age_index(f, r, a);
                pred_comp[static_cast<size_t>(a)] =
                    (pred_capture_age[fra_idx] + Type(1e-12)) /
                    (pred_total + Type(1e-12 * static_cast<double>(a_count)));
              }

              const Type composition_log_likelihood =
                  ddirichlet_multinomial_linear(
                      obs_comp, pred_comp, theta_comp[static_cast<size_t>(f)],
                      true);
              const Type composition_weight =
                  Type(controls_m.composition_likelihood_weight_m);
              nll_composition -=
                  composition_weight * composition_log_likelihood;
              nll -= composition_weight * composition_log_likelihood;

              if (controls_m.report_observation_predictions_m) {
                const std::string row =
                    "y" + std::to_string(y + 1) + "_s" + std::to_string(s + 1) +
                    "_f" + std::to_string(f + 1) + "_r" + std::to_string(r + 1);
                ctx.report("diag_comp_theta_" + row,
                           theta_comp[static_cast<size_t>(f)]);
                for (int a = 0; a < a_count; ++a) {
                  ctx.report("diag_comp_pred_" + row + "_a" +
                                 std::to_string(a + 1),
                             pred_comp[static_cast<size_t>(a)]);
                }
              }
            }

            const double obs_index = data_m.observed_index_m[fysr_idx];
            if (use_index && obs_index > 0.0) {
              const Type pred_index = index_q[static_cast<size_t>(f)] *
                                          pred_vulnerable_bio[fr_idx] +
                                      Type(1e-12);
              nll_index -= dnorm(log(Type(obs_index)), log(pred_index),
                                 sigma_index[static_cast<size_t>(f)], true);
              nll -= dnorm(log(Type(obs_index)), log(pred_index),
                           sigma_index[static_cast<size_t>(f)], true);
              if (controls_m.report_observation_predictions_m) {
                const std::string row =
                    "y" + std::to_string(y + 1) + "_s" + std::to_string(s + 1) +
                    "_f" + std::to_string(f + 1) + "_r" + std::to_string(r + 1);
                ctx.report("diag_index_pred_" + row, pred_index);
                ctx.report("diag_index_sigma_" + row,
                           sigma_index[static_cast<size_t>(f)]);
              }
            }

            const double obs_retained =
                data_m.observed_retained_biomass_m[fysr_idx];
            if (use_retained && obs_retained > 0.0) {
              const Type pred_retained =
                  pred_retained_bio[fr_idx] + Type(1e-12);
              nll_retained -=
                  dnorm(log(Type(obs_retained)), log(pred_retained),
                        sigma_ret_bio[static_cast<size_t>(f)], true);
              nll -= dnorm(log(Type(obs_retained)), log(pred_retained),
                           sigma_ret_bio[static_cast<size_t>(f)], true);
              if (controls_m.report_observation_predictions_m) {
                const std::string row =
                    "y" + std::to_string(y + 1) + "_s" + std::to_string(s + 1) +
                    "_f" + std::to_string(f + 1) + "_r" + std::to_string(r + 1);
                ctx.report("diag_retained_pred_" + row, pred_retained);
                ctx.report("diag_retained_sigma_" + row,
                           sigma_ret_bio[static_cast<size_t>(f)]);
              }
            }

            const double obs_discard =
                data_m.observed_discard_biomass_m[fysr_idx];
            if (use_discard && obs_discard > 0.0) {
              const Type pred_discard = pred_discard_bio[fr_idx] + Type(1e-12);
              nll_discard -= dnorm(log(Type(obs_discard)), log(pred_discard),
                                   sigma_dis_bio[static_cast<size_t>(f)], true);
              nll -= dnorm(log(Type(obs_discard)), log(pred_discard),
                           sigma_dis_bio[static_cast<size_t>(f)], true);
              if (controls_m.report_observation_predictions_m) {
                const std::string row =
                    "y" + std::to_string(y + 1) + "_s" + std::to_string(s + 1) +
                    "_f" + std::to_string(f + 1) + "_r" + std::to_string(r + 1);
                ctx.report("diag_discard_pred_" + row, pred_discard);
                ctx.report("diag_discard_sigma_" + row,
                           sigma_dis_bio[static_cast<size_t>(f)]);
              }
            }
          }
        }

        for (int from = 0; from < r_count; ++from) {
          for (int to = 0; to < r_count; ++to) {
            const Type p = move_probs[movement_probability_index(s, from, to)];

            for (int a = 0; a < a_count; ++a) {
              const size_t from_idx = data_m.region_age_index(from, a);
              const size_t to_idx = data_m.region_age_index(to, a);
              moved[to_idx] += survivors[from_idx] * p;
            }
          }
        }

        state.swap(moved);
      }

      nll_recruitment -= dnorm(recruit_dev[static_cast<size_t>(y)], Type(0.0),
                               sigma_recruit, true);
      nll -= dnorm(recruit_dev[static_cast<size_t>(y)], Type(0.0),
                   sigma_recruit, true);

      if (y + 1 < y_count) {
        const Type ssb_y = spawning_biomass(state, y);
        const Type expected_rec =
            beverton_holt_recruitment(ssb_y, r0, steepness, ssb0);

        std::vector<Type> next_state(state.size(), Type(0.0));

        for (int r = 0; r < r_count; ++r) {
          const Type recruit_share = Type(
              data_m.regional_recruit_proportions_m[static_cast<size_t>(r)]);
          next_state[data_m.region_age_index(r, 0)] =
              expected_rec * recruit_share *
              exp(recruit_dev[static_cast<size_t>(y + 1)] -
                  Type(0.5) * sigma_recruit * sigma_recruit);

          for (int a = 0; a < a_count - 1; ++a) {
            const size_t from_idx = data_m.region_age_index(r, a);
            const size_t to_idx = data_m.region_age_index(r, a + 1);
            next_state[to_idx] += state[from_idx];
          }

          const int plus_age = a_count - 1;
          const size_t plus_idx = data_m.region_age_index(r, plus_age);
          next_state[plus_idx] += state[plus_idx];
        }

        state.swap(next_state);
      }
    }

    ctx.report("r0", r0);
    ctx.report("steepness", steepness);
    ctx.report("sigma_recruit", sigma_recruit);
    ctx.report("nll_prior", nll_prior);
    ctx.report("nll_availability_penalty", nll_availability_penalty);
    ctx.report("nll_movement_penalty", nll_movement_penalty);
    ctx.report("nll_index", nll_index);
    ctx.report("nll_composition", nll_composition);
    ctx.report("nll_retained", nll_retained);
    ctx.report("nll_discard", nll_discard);
    ctx.report("nll_recruitment", nll_recruitment);
    ctx.report("nll_decomposition_total", nll_prior + nll_availability_penalty +
                                              nll_movement_penalty + nll_index +
                                              nll_composition + nll_retained +
                                              nll_discard + nll_recruitment);
    for (int f = 0; f < f_count; ++f) {
      ctx.report("diag_log_q_f" + std::to_string(f + 1),
                 log_q[static_cast<size_t>(f)]);
      ctx.report("diag_q_f" + std::to_string(f + 1), q[static_cast<size_t>(f)]);
      ctx.report("diag_log_index_q_f" + std::to_string(f + 1),
                 log_index_q[static_cast<size_t>(f)]);
      ctx.report("diag_index_q_f" + std::to_string(f + 1),
                 index_q[static_cast<size_t>(f)]);
      for (int s = 0; s < s_count; ++s) {
        for (int r = 0; r < r_count; ++r) {
          const Type log_avail =
              log_availability_scale[fleet_season_region_index(f, s, r)];
          ctx.report("diag_log_availability_f" + std::to_string(f + 1) + "_s" +
                         std::to_string(s + 1) + "_r" + std::to_string(r + 1),
                     log_avail);
          ctx.report("diag_log_q_availability_f" + std::to_string(f + 1) +
                         "_s" + std::to_string(s + 1) + "_r" +
                         std::to_string(r + 1),
                     log_q[static_cast<size_t>(f)] + log_avail);
        }
      }
    }
    ctx.adreport("ssb0", ssb0);
    const Type ssb_terminal = spawning_biomass(state, data_m.n_years_m - 1);
    ctx.adreport("ssb_terminal", ssb_terminal);
    ctx.adreport("depletion_terminal", ssb_terminal / (ssb0 + Type(1e-12)));

    for (int y = 0; y < y_count; ++y) {
      ctx.adreport("ssb_year_" + std::to_string(y + 1),
                   ssb_by_year[static_cast<size_t>(y)]);
    }

    for (int s = 0; s < data_m.n_seasons_m; ++s) {
      for (int from = 0; from < data_m.n_regions_m; ++from) {
        for (int to = 0; to < data_m.n_regions_m; ++to) {
          ctx.adreport("move_prob_s" + std::to_string(s + 1) + "_from_" +
                           std::to_string(from + 1) + "_to_" +
                           std::to_string(to + 1),
                       move_probs[movement_probability_index(s, from, to)]);
        }
      }
    }

    for (int r = 0; r < data_m.n_regions_m; ++r) {
      Type region_bio = Type(0.0);
      for (int a = 0; a < data_m.n_ages_m; ++a) {
        const size_t idx = data_m.region_age_index(r, a);
        ctx.adreport("terminal_numbers_region_" + std::to_string(r + 1) +
                         "_age_" + std::to_string(a + 1),
                     state[idx]);
        const Type w = Type(data_m.weight_at_age_m[data_m.year_age_index(
            data_m.n_years_m - 1, a)]);
        region_bio += state[idx] * w;
      }
      ctx.adreport("terminal_biomass_region_" + std::to_string(r + 1),
                   region_bio);
    }

    return nll;
  }

  TunaSpatialAssessmentData data_m;
  TunaAssessmentControls controls_m;

  int per_fleet_parameter_count() const {
    return 10 + (controls_m.estimate_availability_scales_m
                     ? (controls_m.availability_by_fleet_only_m
                            ? 1
                            : data_m.n_seasons_m * data_m.n_regions_m)
                     : 0);
  }

  int movement_logit_count() const {
    return movement_parameter_season_count() * data_m.n_regions_m *
           (data_m.n_regions_m - 1);
  }

  int movement_parameter_season_count() const {
    return controls_m.share_movement_across_seasons_m ? 1 : data_m.n_seasons_m;
  }

  int movement_parameter_season_index(int season) const {
    return controls_m.share_movement_across_seasons_m ? 0 : season;
  }

  size_t fleet_age_index(int fleet, int age) const {
    return static_cast<size_t>(fleet) * static_cast<size_t>(data_m.n_ages_m) +
           static_cast<size_t>(age);
  }

  size_t fleet_region_index(int fleet, int region) const {
    return static_cast<size_t>(fleet) *
               static_cast<size_t>(data_m.n_regions_m) +
           static_cast<size_t>(region);
  }

  size_t fleet_region_age_index(int fleet, int region, int age) const {
    return (static_cast<size_t>(fleet) *
                static_cast<size_t>(data_m.n_regions_m) +
            static_cast<size_t>(region)) *
               static_cast<size_t>(data_m.n_ages_m) +
           static_cast<size_t>(age);
  }

  size_t fleet_season_region_index(int fleet, int season, int region) const {
    return (static_cast<size_t>(fleet) *
                static_cast<size_t>(data_m.n_seasons_m) +
            static_cast<size_t>(season)) *
               static_cast<size_t>(data_m.n_regions_m) +
           static_cast<size_t>(region);
  }

  size_t movement_logit_index(int season, int from, int to) const {
    return (static_cast<size_t>(season) *
                static_cast<size_t>(data_m.n_regions_m) +
            static_cast<size_t>(from)) *
               static_cast<size_t>(data_m.n_regions_m - 1) +
           static_cast<size_t>(to);
  }

  size_t movement_probability_index(int season, int from, int to) const {
    return (static_cast<size_t>(season) *
                static_cast<size_t>(data_m.n_regions_m) +
            static_cast<size_t>(from)) *
               static_cast<size_t>(data_m.n_regions_m) +
           static_cast<size_t>(to);
  }

  double movement_initial_logit(int season, int from, int to) const {
    const int baseline_to = data_m.n_regions_m - 1;
    const double p_to =
        data_m.movement_matrix_m[data_m.season_region_region_index(season, from,
                                                                   to)];
    const double p_base =
        data_m.movement_matrix_m[data_m.season_region_region_index(
            season, from, baseline_to)];
    return std::log((p_to + 1e-12) / (p_base + 1e-12));
  }

  double movement_initial_logit_shared(int from, int to) const {
    const int baseline_to = data_m.n_regions_m - 1;

    double p_to_sum = 0.0;
    double p_base_sum = 0.0;
    for (int s = 0; s < data_m.n_seasons_m; ++s) {
      p_to_sum += data_m.movement_matrix_m[data_m.season_region_region_index(
          s, from, to)];
      p_base_sum += data_m.movement_matrix_m[data_m.season_region_region_index(
          s, from, baseline_to)];
    }

    const double p_to = p_to_sum / static_cast<double>(data_m.n_seasons_m);
    const double p_base = p_base_sum / static_cast<double>(data_m.n_seasons_m);
    return std::log((p_to + 1e-12) / (p_base + 1e-12));
  }

  template <typename Type>
  Type effective_availability(int fleet, int season, int region, int age,
                              const Type &availability_scale) const {
    Type surface = Type(1.0);
    if (!data_m.availability_surface_m.empty()) {
      surface = Type(
          data_m.availability_surface_m[data_m.fleet_season_region_age_index(
              fleet, season, region, age)]);
    }
    return surface * availability_scale;
  }

  template <typename Type> Type bounded_age50(const Type &raw) const {
    const Type unit = Type(1.0) / (Type(1.0) + exp(-raw));
    return Type(1.0) + Type(data_m.n_ages_m - 1) * unit;
  }

  template <typename Type>
  Type logistic_age(const Type &age, const Type &a50, const Type &slope) const {
    const Type scaled = slope * (age - a50);
    return Type(1.0) / (Type(1.0) + exp(-scaled));
  }

  template <typename Type>
  Type normalized_logistic_age(const Type &age, const Type &a50,
                               const Type &slope) const {
    const Type oldest = Type(data_m.n_ages_m);
    return logistic_age(age, a50, slope) /
           (logistic_age(oldest, a50, slope) + Type(1e-12));
  }

  template <typename Type>
  Type beverton_holt_recruitment(const Type &ssb, const Type &r0, const Type &h,
                                 const Type &ssb0) const {
    const Type four_h_r0 = Type(4.0) * h * r0;
    const Type one_minus_h = Type(1.0) - h;

    const Type num = four_h_r0 * ssb;
    const Type den =
        ssb0 * one_minus_h + ssb * (Type(5.0) * h - Type(1.0)) + Type(1e-12);

    return num / den;
  }

  template <typename Type>
  std::vector<Type>
  initialize_state(const Type &r0, const Type &sigma_recruit,
                   const std::vector<Type> &recruit_dev) const {
    std::vector<Type> state(static_cast<size_t>(data_m.n_regions_m) *
                                static_cast<size_t>(data_m.n_ages_m),
                            Type(0.0));

    for (int r = 0; r < data_m.n_regions_m; ++r) {
      const Type recruit_share =
          Type(data_m.regional_recruit_proportions_m[static_cast<size_t>(r)]);
      state[data_m.region_age_index(r, 0)] =
          r0 * recruit_share *
          exp(recruit_dev[0] - Type(0.5) * sigma_recruit * sigma_recruit);

      for (int a = 1; a < data_m.n_ages_m - 1; ++a) {
        const size_t prev_idx = data_m.region_age_index(r, a - 1);
        const size_t idx = data_m.region_age_index(r, a);
        const Type m =
            Type(data_m.natural_mortality_at_age_m[static_cast<size_t>(a - 1)]);
        state[idx] = state[prev_idx] * exp(-m);
      }

      const int plus_age = data_m.n_ages_m - 1;
      const size_t prev_idx = data_m.region_age_index(r, plus_age - 1);
      const size_t plus_idx = data_m.region_age_index(r, plus_age);
      const Type m_plus = Type(
          data_m.natural_mortality_at_age_m[static_cast<size_t>(plus_age)]);
      const Type surv_plus = exp(-m_plus);
      state[plus_idx] =
          state[prev_idx] * surv_plus / (Type(1.0) - surv_plus + Type(1e-12));
    }

    return state;
  }

  template <typename Type>
  Type spawning_biomass(const std::vector<Type> &state, int year) const {
    Type out = Type(0.0);

    for (int r = 0; r < data_m.n_regions_m; ++r) {
      for (int a = 0; a < data_m.n_ages_m; ++a) {
        const size_t ra_idx = data_m.region_age_index(r, a);
        const Type n = state[ra_idx];
        const Type w =
            Type(data_m.weight_at_age_m[data_m.year_age_index(year, a)]);
        const Type mat = Type(data_m.maturity_at_age_m[static_cast<size_t>(a)]);
        const Type m =
            Type(data_m.natural_mortality_at_age_m[static_cast<size_t>(a)]);
        const Type surv = exp(-m * Type(data_m.spawning_fraction_m));

        out += n * w * mat * surv;
      }
    }

    return out;
  }
};

} // namespace quadra

#endif // QUADRA_TUNA_SPATIAL_ASSESSMENT_MODEL_HPP
