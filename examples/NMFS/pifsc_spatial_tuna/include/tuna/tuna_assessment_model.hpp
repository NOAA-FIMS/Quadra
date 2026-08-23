#ifndef QUADRA_TUNA_ASSESSMENT_MODEL_HPP
#define QUADRA_TUNA_ASSESSMENT_MODEL_HPP
#pragma once

#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "core/model/parameter.hpp"
#include "core/model/quadra_model.hpp"
#include "data.hpp"
#include "math/distributions.hpp"

namespace quadra {

class AdvancedTunaAssessmentModel
    : public QuadraModel<AdvancedTunaAssessmentModel> {
public:
  explicit AdvancedTunaAssessmentModel(TunaAssessmentData data)
      : data_m(std::move(data)) {
    data_m.validate();
  }

  ParameterSet parameter_set() const {
    ParameterSet p;
    p.add("log_r0", 10.0);
    p.add("logit_steepness", 0.0);

    for (int f = 0; f < data_m.n_fleets_m; ++f) {
      p.add("log_q_fleet_" + std::to_string(f + 1), -8.0);
      p.add("log_sel50_fleet_" + std::to_string(f + 1),
            std::log(0.5 * static_cast<double>(data_m.n_ages_m)));
      p.add("log_sel_slope_fleet_" + std::to_string(f + 1), 0.0);
      p.add("log_sigma_index_fleet_" + std::to_string(f + 1), -1.0);
      p.add("log_theta_comp_fleet_" + std::to_string(f + 1), -2.0);
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
    const size_t n_fixed = static_cast<size_t>(2 + 5 * data_m.n_fleets_m + 1);
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
    const int y_count = data_m.n_years_m;
    const int a_count = data_m.n_ages_m;
    const int f_count = data_m.n_fleets_m;

    const size_t expected = parameter_set().size();
    if (parameters.size() != expected) {
      throw std::invalid_argument(
          "AdvancedTunaAssessmentModel: parameter vector length mismatch");
    }

    size_t pos = 0;
    const Type log_r0 = parameters[pos++];
    const Type logit_h = parameters[pos++];

    std::vector<Type> log_q(static_cast<size_t>(f_count));
    std::vector<Type> log_sel50(static_cast<size_t>(f_count));
    std::vector<Type> log_sel_slope(static_cast<size_t>(f_count));
    std::vector<Type> log_sigma_index(static_cast<size_t>(f_count));
    std::vector<Type> log_theta_comp(static_cast<size_t>(f_count));

    for (int f = 0; f < f_count; ++f) {
      const size_t ff = static_cast<size_t>(f);
      log_q[ff] = parameters[pos++];
      log_sel50[ff] = parameters[pos++];
      log_sel_slope[ff] = parameters[pos++];
      log_sigma_index[ff] = parameters[pos++];
      log_theta_comp[ff] = parameters[pos++];
    }

    const Type log_sigma_recruit = parameters[pos++];

    std::vector<Type> recruit_dev(static_cast<size_t>(y_count));
    for (int y = 0; y < y_count; ++y) {
      recruit_dev[static_cast<size_t>(y)] = parameters[pos++];
    }

    Type nll = Type(0.0);

    const Type r0 = exp(log_r0);
    const Type steepness = Type(0.2) + Type(0.8) * logistic(logit_h);
    const Type sigma_recruit = exp(log_sigma_recruit);

    std::vector<Type> q(static_cast<size_t>(f_count));
    std::vector<Type> sel50(static_cast<size_t>(f_count));
    std::vector<Type> sel_slope(static_cast<size_t>(f_count));
    std::vector<Type> sigma_index(static_cast<size_t>(f_count));
    std::vector<Type> theta_comp(static_cast<size_t>(f_count));

    for (int f = 0; f < f_count; ++f) {
      const size_t ff = static_cast<size_t>(f);
      q[ff] = exp(log_q[ff]);
      sel50[ff] = exp(log_sel50[ff]);
      sel_slope[ff] = exp(log_sel_slope[ff]);
      sigma_index[ff] = exp(log_sigma_index[ff]);
      theta_comp[ff] = exp(log_theta_comp[ff]);
    }

    const size_t n_ya =
        static_cast<size_t>(y_count) * static_cast<size_t>(a_count);
    std::vector<Type> numbers_at_age(n_ya, Type(0.0));

    numbers_at_age[data_m.year_age_index(0, 0)] =
        r0 * exp(recruit_dev[0] - Type(0.5) * sigma_recruit * sigma_recruit);

    for (int a = 1; a < a_count - 1; ++a) {
      const size_t prev_idx = data_m.year_age_index(0, a - 1);
      const size_t idx = data_m.year_age_index(0, a);
      const Type m_prev =
          Type(data_m.natural_mortality_at_age_m[static_cast<size_t>(a - 1)]);
      numbers_at_age[idx] = numbers_at_age[prev_idx] * exp(-m_prev);
    }

    {
      const int plus_age = a_count - 1;
      const size_t prev_idx = data_m.year_age_index(0, plus_age - 1);
      const size_t plus_idx = data_m.year_age_index(0, plus_age);
      const Type m_plus = Type(
          data_m.natural_mortality_at_age_m[static_cast<size_t>(plus_age)]);
      const Type surv_plus = exp(-m_plus);
      const Type denom = Type(1.0) - surv_plus + Type(1e-12);
      numbers_at_age[plus_idx] = numbers_at_age[prev_idx] * surv_plus / denom;
    }

    const Type ssb0 = spawning_biomass(numbers_at_age, 0);

    std::vector<Type> selectivity(static_cast<size_t>(f_count) *
                                  static_cast<size_t>(a_count));
    for (int f = 0; f < f_count; ++f) {
      for (int a = 0; a < a_count; ++a) {
        const size_t idx =
            static_cast<size_t>(f) * static_cast<size_t>(a_count) +
            static_cast<size_t>(a);
        selectivity[idx] = logistic_age(Type(static_cast<double>(a) + 1.0),
                                        sel50[static_cast<size_t>(f)],
                                        sel_slope[static_cast<size_t>(f)]);
      }
    }

    for (int y = 0; y < y_count; ++y) {
      for (int f = 0; f < f_count; ++f) {
        const size_t fy_idx = data_m.fleet_year_index(f, y);
        const double obs_index = data_m.observed_index_m[fy_idx];

        std::vector<Type> pred_catch_age(static_cast<size_t>(a_count),
                                         Type(0.0));
        std::vector<Type> pred_catch_prop(static_cast<size_t>(a_count),
                                          Type(0.0));
        Type pred_catch_sum = Type(0.0);
        Type pred_vulnerable_bio = Type(0.0);

        for (int a = 0; a < a_count; ++a) {
          const size_t ya_idx = data_m.year_age_index(y, a);
          const Type n = numbers_at_age[ya_idx];
          const Type m =
              Type(data_m.natural_mortality_at_age_m[static_cast<size_t>(a)]);

          Type f_total = Type(0.0);
          for (int ff = 0; ff < f_count; ++ff) {
            const size_t ffy_idx = data_m.fleet_year_index(ff, y);
            const Type effort = Type(data_m.effort_m[ffy_idx]);
            const size_t sel_idx =
                static_cast<size_t>(ff) * static_cast<size_t>(a_count) +
                static_cast<size_t>(a);
            f_total +=
                q[static_cast<size_t>(ff)] * effort * selectivity[sel_idx];
          }

          const size_t sel_idx =
              static_cast<size_t>(f) * static_cast<size_t>(a_count) +
              static_cast<size_t>(a);
          const Type effort_f = Type(data_m.effort_m[fy_idx]);
          const Type f_fleet =
              q[static_cast<size_t>(f)] * effort_f * selectivity[sel_idx];

          const Type z = m + f_total + Type(1e-12);
          const Type catch_n = n * (f_fleet / z) * (Type(1.0) - exp(-z));
          pred_catch_age[static_cast<size_t>(a)] = catch_n;
          pred_catch_sum += catch_n;

          const Type w = Type(data_m.weight_at_age_m[ya_idx]);
          pred_vulnerable_bio += n * w * selectivity[sel_idx];
        }

        const std::vector<int> obs_comp = observed_catch_for_fleet_year(f, y);
        int obs_total = 0;
        for (int a = 0; a < a_count; ++a) {
          obs_total += obs_comp[static_cast<size_t>(a)];
        }

        if (obs_total > 0) {
          for (int a = 0; a < a_count; ++a) {
            pred_catch_prop[static_cast<size_t>(a)] =
                (pred_catch_age[static_cast<size_t>(a)] + Type(1e-12)) /
                (pred_catch_sum + Type(1e-12 * static_cast<double>(a_count)));
          }

          nll -= ddirichlet_multinomial_linear(
              obs_comp, pred_catch_prop, theta_comp[static_cast<size_t>(f)],
              true);
        }

        if (obs_index > 0.0) {
          const Type pred_index =
              q[static_cast<size_t>(f)] * pred_vulnerable_bio + Type(1e-12);

          nll -= dnorm(log(Type(obs_index)), log(pred_index),
                       sigma_index[static_cast<size_t>(f)], true);
        }
      }

      nll -= dnorm(recruit_dev[static_cast<size_t>(y)], Type(0.0),
                   sigma_recruit, true);

      if (y + 1 < y_count) {
        const Type ssb_y = spawning_biomass(numbers_at_age, y);
        const Type expected_rec =
            beverton_holt_recruitment(ssb_y, r0, steepness, ssb0);

        const size_t rec_idx = data_m.year_age_index(y + 1, 0);
        numbers_at_age[rec_idx] =
            expected_rec * exp(recruit_dev[static_cast<size_t>(y + 1)] -
                               Type(0.5) * sigma_recruit * sigma_recruit);

        for (int a = 0; a < a_count - 1; ++a) {
          const size_t current_idx = data_m.year_age_index(y, a);
          const size_t next_idx = data_m.year_age_index(y + 1, a + 1);

          Type f_total = Type(0.0);
          for (int ff = 0; ff < f_count; ++ff) {
            const size_t ffy_idx = data_m.fleet_year_index(ff, y);
            const Type effort = Type(data_m.effort_m[ffy_idx]);
            const size_t sel_idx =
                static_cast<size_t>(ff) * static_cast<size_t>(a_count) +
                static_cast<size_t>(a);
            f_total +=
                q[static_cast<size_t>(ff)] * effort * selectivity[sel_idx];
          }

          const Type m =
              Type(data_m.natural_mortality_at_age_m[static_cast<size_t>(a)]);
          const Type z = m + f_total;
          numbers_at_age[next_idx] = numbers_at_age[current_idx] * exp(-z);
        }

        {
          const int plus_age = a_count - 1;
          const size_t plus_from = data_m.year_age_index(y, plus_age);
          const size_t plus_to = data_m.year_age_index(y + 1, plus_age);

          Type f_total = Type(0.0);
          for (int ff = 0; ff < f_count; ++ff) {
            const size_t ffy_idx = data_m.fleet_year_index(ff, y);
            const Type effort = Type(data_m.effort_m[ffy_idx]);
            const size_t sel_idx =
                static_cast<size_t>(ff) * static_cast<size_t>(a_count) +
                static_cast<size_t>(plus_age);
            f_total +=
                q[static_cast<size_t>(ff)] * effort * selectivity[sel_idx];
          }

          const Type m_plus = Type(
              data_m.natural_mortality_at_age_m[static_cast<size_t>(plus_age)]);
          const Type z_plus = m_plus + f_total;
          numbers_at_age[plus_to] += numbers_at_age[plus_from] * exp(-z_plus);
        }
      }
    }

    ctx.report("r0", r0);
    ctx.report("steepness", steepness);
    ctx.report("sigma_recruit", sigma_recruit);
    ctx.adreport("ssb0", ssb0);
    ctx.adreport("ssb_terminal",
                 spawning_biomass(numbers_at_age, data_m.n_years_m - 1));

    return nll;
  }

private:
  TunaAssessmentData data_m;

  template <typename Type>
  Type logistic_age(const Type &age, const Type &a50, const Type &slope) const {
    const Type scaled = (age - a50) / (slope + Type(1e-12));
    return Type(1.0) / (Type(1.0) + exp(-scaled));
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
  Type spawning_biomass(const std::vector<Type> &numbers_at_age,
                        int year) const {
    Type out = Type(0.0);

    for (int a = 0; a < data_m.n_ages_m; ++a) {
      const size_t ya_idx = data_m.year_age_index(year, a);
      const Type n = numbers_at_age[ya_idx];
      const Type w = Type(data_m.weight_at_age_m[ya_idx]);
      const Type mat = Type(data_m.maturity_at_age_m[static_cast<size_t>(a)]);
      const Type m =
          Type(data_m.natural_mortality_at_age_m[static_cast<size_t>(a)]);
      const Type surv = exp(-m * Type(data_m.spawning_fraction_m));

      out += n * w * mat * surv;
    }

    return out;
  }

  std::vector<int> observed_catch_for_fleet_year(int fleet, int year) const {
    std::vector<int> out(static_cast<size_t>(data_m.n_ages_m), 0);

    for (int a = 0; a < data_m.n_ages_m; ++a) {
      const size_t idx = data_m.fleet_year_age_index(fleet, year, a);
      out[static_cast<size_t>(a)] = data_m.observed_catch_numbers_m[idx];
    }

    return out;
  }
};

} // namespace quadra

#endif // QUADRA_TUNA_ASSESSMENT_MODEL_HPP
