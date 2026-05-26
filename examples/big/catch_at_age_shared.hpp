#pragma once

// examples/big/catch_at_age_laplace.cpp
//
// Big fisheries-oriented Quadra example:
//   age-structured catch-at-age model with recruitment random effects.
//
// This is an example/benchmark model, not a production assessment model.

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

#include "../../core/laplace/laplace_evaluator.hpp"
#include "../../core/model/parameter.hpp"

#include "../../core/optimizer.hpp"

// quadra_big_example_plateau_cap_v2
// Trace evidence: objective is effectively flat by ~150 outer evaluations.
// This keeps the example from grinding after the fitted solution is reached.
#ifndef QUADRA_BIG_EXAMPLE_MAX_ITER
#define QUADRA_BIG_EXAMPLE_MAX_ITER 150
#endif
#ifndef QUADRA_BIG_EXAMPLE_EPSILON
#define QUADRA_BIG_EXAMPLE_EPSILON 1e-4
#endif

DECLARE_ADGRAPH();

namespace example {

// quadra_big_example_derived_reporting_v1
template <typename Result, typename = void>
struct has_final_random_effects : std::false_type {};

template <typename Result>
struct has_final_random_effects<
    Result, std::void_t<decltype(std::declval<Result>().final_random_effects)>>
    : std::true_type {};

template <typename T> T square(const T &x) { return x * x; }

template <typename T> T inv_logit(const T &x) {
  return T(1.0) / (T(1.0) + exp(-x));
}

template <typename T>
T logistic_selectivity(const T &age, const T &a50, const T &slope) {
  return T(1.0) / (T(1.0) + exp(-slope * (age - a50)));
}

template <typename T> T safe_log(const T &x) { return log(x + T(1.0e-12)); }

// quadra_big_example_stabilized_objective

template <typename T>
T normal_prior_nll(const T &x, const double mean, const double sd) {
  const T z = (x - T(mean)) / T(sd);
  return T(0.5) * z * z;
}

template <typename T> T positive_floor(const T &x, const double floor) {
  return T(floor) + exp(x);
}

struct CatchAtAgeData {
  int n_years = 30;
  int n_ages = 8;

  std::vector<double> index_obs;
  std::vector<double> catch_obs;
  std::vector<std::vector<double>> age_comp_obs;

  CatchAtAgeData() {
    index_obs.resize(static_cast<std::size_t>(n_years));
    catch_obs.resize(static_cast<std::size_t>(n_years));
    age_comp_obs.assign(
        static_cast<std::size_t>(n_years),
        std::vector<double>(static_cast<std::size_t>(n_ages), 0.0));

    for (int y = 0; y < n_years; ++y) {
      const double trend = std::exp(-0.025 * static_cast<double>(y));
      const double pulse = 1.0 + 0.12 * std::sin(0.45 * static_cast<double>(y));

      index_obs[static_cast<std::size_t>(y)] = 850.0 * trend * pulse + 20.0;
      catch_obs[static_cast<std::size_t>(y)] =
          120.0 * pulse + 5.0 * std::cos(0.2 * y);

      double total = 0.0;
      for (int a = 0; a < n_ages; ++a) {
        const double age = static_cast<double>(a + 1);
        const double z = (age - 4.5 - 0.015 * y) / 1.65;
        const double val = std::exp(-0.5 * z * z) + 0.02;
        age_comp_obs[static_cast<std::size_t>(y)][static_cast<std::size_t>(a)] =
            val;
        total += val;
      }

      for (int a = 0; a < n_ages; ++a) {
        age_comp_obs[static_cast<std::size_t>(y)]
                    [static_cast<std::size_t>(a)] /= total;
      }
    }
  }
};

struct CatchAtAgeLaplaceModel {
  CatchAtAgeData data;

  template <typename Context> void initialize(Context &) {
    // This example does not emit reports, but Quadra's model interface
    // expects an initialize(ctx) hook.
  }

  template <typename T, typename Context>
  T evaluate(const std::vector<T> &p, Context &) const {
    return evaluate_impl<T>(p);
  }

  template <typename T> T evaluate(const std::vector<T> &p) const {
    return evaluate_impl<T>(p);
  }

  template <typename U, typename V>
  void trace_add_nll(U &nll, const V &term, const char *label) const {
    nll += term;
#ifdef QUADRA_TRACE_MODEL_OBJECTIVE
    if constexpr (std::is_same_v<U, double>) {
      std::cout << "[objective trace] " << label << " term=" << term
                << " cumulative=" << nll << std::endl;
    }
#endif
  }

  template <typename T> T evaluate_impl(const std::vector<T> &p) const {
    const int n_years = data.n_years;
    const int n_ages = data.n_ages;

    const T log_R0 = p[0];
    const T log_M = p[1];
    const T log_q = p[2];
    const T log_Fbar = p[3];
    const T sel50_raw = p[4];
    const T log_sel_slope = p[5];
    const T log_sigma_index = p[6];
    const T log_sigma_catch = p[7];
    const T log_sigma_rec = p[8];

    const T R0 = exp(log_R0);
    const T M = exp(log_M);
    const T q = exp(log_q);
    const T Fbar = exp(log_Fbar);
    const T sel50 = T(1.0) + T(n_ages - 1) * inv_logit(sel50_raw);
    const T sel_slope = exp(log_sel_slope);
    const T sigma_index = positive_floor(log_sigma_index, 0.03);
    const T sigma_catch = positive_floor(
        log_sigma_catch, 0.08); // quadra_big_example_model_improvement_v2
    const T sigma_rec = positive_floor(log_sigma_rec, 0.05);

    std::vector<T> N(static_cast<std::size_t>(n_ages), T(0.0));
    std::vector<T> N_next(static_cast<std::size_t>(n_ages), T(0.0));
    std::vector<T> sel(static_cast<std::size_t>(n_ages), T(0.0));
    std::vector<T> catch_at_age(static_cast<std::size_t>(n_ages), T(0.0));

    for (int a = 0; a < n_ages; ++a) {
      const T age = T(a + 1);
      sel[static_cast<std::size_t>(a)] =
          logistic_selectivity(age, sel50, sel_slope);
      N[static_cast<std::size_t>(a)] = R0 * exp(-M * T(a));
    }

    T nll = T(0.0);

    // quadra_big_example_model_improvement_v1
    // Recruitment deviations and log_R0 can trade off. A weak
    // sum-to-zero-style penalty improves identifiability and makes the
    // demo optimization surface less sloppy.
    T mean_rec_dev = T(0.0);
    for (int yy = 0; yy < n_years; ++yy) {
      mean_rec_dev += p[9 + yy];
    }
    mean_rec_dev /= T(n_years);
    trace_add_nll(nll, T(0.5) * square(mean_rec_dev / T(0.10)),
                  "recruitment_mean_penalty");

    // v2: tighter priors reduce R0-M-q-Fbar confounding in the demo.
    trace_add_nll(nll, normal_prior_nll(log_R0, std::log(900.0), 0.55),
                  "prior_log_R0");
    trace_add_nll(nll, normal_prior_nll(log_M, std::log(0.25), 0.22),
                  "prior_log_M");
    trace_add_nll(nll, normal_prior_nll(log_q, std::log(0.15), 0.55),
                  "prior_log_q");
    trace_add_nll(nll, normal_prior_nll(log_Fbar, std::log(0.18), 0.40),
                  "prior_log_Fbar");
    trace_add_nll(nll, normal_prior_nll(sel50_raw, 0.0, 1.50),
                  "prior_sel50_raw");
    trace_add_nll(nll, normal_prior_nll(log_sel_slope, std::log(1.25), 0.75),
                  "prior_log_sel_slope");

    trace_add_nll(nll, normal_prior_nll(log_sigma_index, std::log(0.20), 0.75),
                  "prior_log_sigma_index");
    trace_add_nll(nll, normal_prior_nll(log_sigma_catch, std::log(0.18), 0.35),
                  "prior_log_sigma_catch");
    trace_add_nll(nll, normal_prior_nll(log_sigma_rec, std::log(0.35), 0.75),
                  "prior_log_sigma_rec");

    for (int y = 0; y < n_years; ++y) {
      const T rec_dev = p[9 + y];

      trace_add_nll(nll, T(0.5) * square(rec_dev / sigma_rec) + log(sigma_rec),
                    "recruitment_prior");

      // Soft stabilizer against pathological random-effect modes.
      trace_add_nll(nll, T(0.5) * square(rec_dev / T(4.0)),
                    "recruitment_soft_stabilizer");

      const T recruitment = R0 * exp(rec_dev - T(0.5) * square(sigma_rec));

      T vulnerable_biomass = T(0.0);
      T total_catch = T(0.0);

      for (int a = 0; a < n_ages; ++a) {
        const T age = T(a + 1);
        const T weight = (age * age * age) / T(100.0);
        const T F_a = Fbar * sel[static_cast<std::size_t>(a)];
        const T Z_a = M + F_a;
        const T C_a =
            N[static_cast<std::size_t>(a)] * (F_a / Z_a) * (T(1.0) - exp(-Z_a));

        catch_at_age[static_cast<std::size_t>(a)] = C_a;
        total_catch += C_a * weight;
        vulnerable_biomass += N[static_cast<std::size_t>(a)] * weight *
                              sel[static_cast<std::size_t>(a)];
      }

      const T pred_index = q * vulnerable_biomass + T(1.0e-9);
      const T pred_catch = total_catch + T(1.0e-9);
      const T obs_index = T(data.index_obs[static_cast<std::size_t>(y)]);
      const T obs_catch = T(data.catch_obs[static_cast<std::size_t>(y)]);

      nll += T(0.5) * square((safe_log(obs_index) - safe_log(pred_index)) /
                             sigma_index) +
             log_sigma_index;
      nll += T(0.5) * square((safe_log(obs_catch) - safe_log(pred_catch)) /
                             sigma_catch) +
             log_sigma_catch;

      const T catch_sum =
          std::accumulate(catch_at_age.begin(), catch_at_age.end(), T(0.0)) +
          T(1.0e-12);

      const T comp_weight = T(80.0);

      std::fill(N_next.begin(), N_next.end(), T(0.0));
      N_next[0] = recruitment;

      for (int a = 1; a < n_ages; ++a) {
        const T F_prev = Fbar * sel[static_cast<std::size_t>(a - 1)];
        const T Z_prev = M + F_prev;
        N_next[static_cast<std::size_t>(a)] =
            N[static_cast<std::size_t>(a - 1)] * exp(-Z_prev);
      }

      const T F_plus = Fbar * sel[static_cast<std::size_t>(n_ages - 1)];
      const T Z_plus = M + F_plus;
      N_next[static_cast<std::size_t>(n_ages - 1)] +=
          N[static_cast<std::size_t>(n_ages - 1)] * exp(-Z_plus);

      N.swap(N_next);
    }

    // Weak regularization for optimizer stability.

    return nll;
  }

  template <typename T> T operator()(const std::vector<T> &p) const {
    return evaluate_impl<T>(p);
  }
};

template <typename T, typename = void>
struct has_objective_value_m : std::false_type {};

template <typename T>
struct has_objective_value_m<
    T, std::void_t<decltype(std::declval<T>().objective_value_m)>>
    : std::true_type {};

template <typename T, typename = void> struct has_value_m : std::false_type {};

template <typename T>
struct has_value_m<T, std::void_t<decltype(std::declval<T>().value_m)>>
    : std::true_type {};

template <typename T, typename = void> struct has_value : std::false_type {};

template <typename T>
struct has_value<T, std::void_t<decltype(std::declval<T>().value)>>
    : std::true_type {};

template <typename Result> double objective_from_result(const Result &result) {
  if constexpr (has_objective_value_m<Result>::value) {
    return result.objective_value_m;
  } else if constexpr (has_value_m<Result>::value) {
    return result.value_m;
  } else if constexpr (has_value<Result>::value) {
    return result.value;
  } else {
    return 0.0;
  }
}

void print_fixed_parameter_report(const std::vector<double> &theta) {
  if (theta.size() < 9) {
    std::cout
        << "fixed parameter report unavailable: expected 9 fixed effects, got "
        << theta.size() << "\n";
    return;
  }

  const double log_R0 = theta[0];
  const double log_M = theta[1];
  const double log_q = theta[2];
  const double log_Fbar = theta[3];
  const double sel50_raw = theta[4];
  const double log_sel_slope = theta[5];
  const double log_sigma_index = theta[6];
  const double log_sigma_catch = theta[7];
  const double log_sigma_rec = theta[8];

  std::cout << "Fixed - effect report" << std::endl;

  std::cout << "R0: " << std::exp(log_R0) << std::endl;
  std::cout << "M: " << std::exp(log_M) << std::endl;

  std::cout << "q: " << std::exp(log_q) << std::endl;
  std::cout << "Fbar: " << std::exp(log_Fbar) << std::endl;
  std::cout << "sel50_raw: " << sel50_raw << std::endl;

  std::cout << "sel_slope: " << std::exp(log_sel_slope) << std::endl;
  std::cout << "sigma_index_floor_adjusted: "
            << 0.03 + std::exp(log_sigma_index) << std::endl;

  const double sigma_catch_report = 0.08 + std::exp(log_sigma_catch);

  std::cout << "sigma_catch_floor_adjusted: " << sigma_catch_report << "\n";

  // quadra_sync_sigma_catch_report_v2

  std::cout << "sigma_catch_excess_over_floor: " << sigma_catch_report - 0.08
            << "\n";

  if (sigma_catch_report - 0.08 < 1.0e-3)

  {

    std::cout << "sigma_catch_note: near floor; catch likelihood is still "
                 "highly informative in this demo.\n";
  }
  std::cout << "sigma_rec_floor_adjusted: " << 0.05 + std::exp(log_sigma_rec)
            << std::endl;
}

template <typename Result>
void print_derived_quantity_report(
    const CatchAtAgeLaplaceModel &model, const std::vector<double> &theta,
    const Result &result, const std::vector<double> &final_random_effects) {
  std::cout << std::endl;
  std::cout << "Derived quantity report" << std::endl;

  if (theta.size() < 9) {
    std::cout << "derived report unavailable: expected 9 fixed effects"
              << std::endl;
    return;
  }
  {
    const auto &u = final_random_effects;

    if (u.size() < static_cast<std::size_t>(model.data.n_years)) {
      std::cout << "derived report unavailable: final_random_effects shorter "
                   "than n_years"
                << std::endl;
      return;
    }

    const int n_years = model.data.n_years;
    const int n_ages = model.data.n_ages;

    const double R0 = std::exp(theta[0]);
    const double M = std::exp(theta[1]);
    const double q = std::exp(theta[2]);
    const double Fbar = std::exp(theta[3]);
    const double sel50_raw = theta[4];
    const double sel_slope = std::exp(theta[5]);

    std::vector<double> selectivity(static_cast<std::size_t>(n_ages), 0.0);
    for (int a = 0; a < n_ages; ++a) {
      const double age = static_cast<double>(a + 1);
      const double sel50 = 3.0 + sel50_raw;
      selectivity[static_cast<std::size_t>(a)] =
          1.0 / (1.0 + std::exp(-(age - sel50) / sel_slope));
    }

    std::vector<double> numbers(static_cast<std::size_t>(n_ages),
                                R0 / static_cast<double>(n_ages));

    std::cout << "Selectivity at age" << std::endl;
    for (int a = 0; a < n_ages; ++a) {
      std::cout << "  age " << std::setw(2) << (a + 1) << ": "
                << selectivity[static_cast<std::size_t>(a)] << std::endl;
    }

    std::cout << std::endl;
    std::cout << "Yearly derived quantities" << std::endl;
    std::cout << std::setw(6) << "year" << std::setw(16) << "recruitment"
              << std::setw(16) << "ssb_proxy" << std::setw(18) << "vuln_biomass"
              << std::setw(16) << "catch_bio" << std::setw(12) << "Fbar"
              << std::setw(14) << "depletion" << std::endl;

    double ssb0 = 0.0;

    for (int y = 0; y < n_years; ++y) {
      const double rec = R0 * std::exp(u[static_cast<std::size_t>(y)]);
      numbers[0] = rec;

      double ssb = 0.0;
      double vulnerable_biomass = 0.0;
      double catch_biomass = 0.0;
      const double Fy = Fbar;

      for (int a = 0; a < n_ages; ++a) {
        const double age = static_cast<double>(a + 1);
        const double weight = (age * age * age) / 100.0;
        const double maturity = 1.0 / (1.0 + std::exp(-(age - 4.0)));

        const double Na = numbers[static_cast<std::size_t>(a)];
        const double Sa = selectivity[static_cast<std::size_t>(a)];
        const double F_at_age = Fy * Sa;
        const double Z = M + F_at_age;

        ssb += Na * weight * maturity;
        vulnerable_biomass += Na * weight * Sa;

        if (Z > 1.0e-12) {
          catch_biomass += Na * weight * (F_at_age / Z) * (1.0 - std::exp(-Z));
        }
      }

      if (y == 0) {
        ssb0 = ssb;
      }

      const double depletion =
          ssb0 > 0.0 ? ssb / ssb0 : std::numeric_limits<double>::quiet_NaN();

      std::cout << std::setw(6) << (y + 1) << std::setw(16) << rec
                << std::setw(16) << ssb << std::setw(18) << vulnerable_biomass
                << std::setw(16) << catch_biomass << std::setw(12) << Fy
                << std::setw(14) << depletion << std::endl;

      for (int a = n_ages - 1; a >= 1; --a) {
        const double Sa_prev = selectivity[static_cast<std::size_t>(a - 1)];
        const double Z_prev = M + Fy * Sa_prev;
        numbers[static_cast<std::size_t>(a)] =
            numbers[static_cast<std::size_t>(a - 1)] * std::exp(-Z_prev);
      }
    }

    std::cout << std::endl;
    std::cout << "Derived report notes" << std::endl;
    std::cout << "  ssb_proxy uses maturity-at-age and cubic weight-at-age "
                 "from the demo model."
              << std::endl;
    std::cout << "  depletion is relative to year-1 ssb_proxy, not an unfished "
                 "equilibrium reference point."
              << std::endl;
    std::cout << "  q is estimated but only affects the index likelihood, not "
                 "the derived biomass scale directly."
              << std::endl;
    (void)q;
  }
}

// quadra_big_example_objective_decomposition_v1
// quadra_fix_objective_decomp_age_comp_obs_v1
struct ObjectiveComponents {
  double fixed_priors = 0.0;
  double recruitment_prior = 0.0;
  double recruitment_mean_penalty = 0.0;
  double recruitment_soft_stabilizer = 0.0;
  double index_likelihood = 0.0;
  double catch_likelihood = 0.0;
  double composition_proxy = 0.0;

  double total() const {
    return fixed_priors + recruitment_prior + recruitment_mean_penalty +
           recruitment_soft_stabilizer + index_likelihood + catch_likelihood +
           composition_proxy;
  }
};

inline double square_double(const double x) { return x * x; }

inline double safe_log_double(const double x) { return std::log(x + 1.0e-12); }

inline double normal_prior_nll_double(const double x, const double mean,
                                      const double sd) {
  const double z = (x - mean) / sd;
  return 0.5 * z * z;
}

inline ObjectiveComponents
evaluate_objective_components(const CatchAtAgeLaplaceModel &model,
                              const std::vector<double> &theta,
                              const std::vector<double> &u) {
  ObjectiveComponents c;

  if (theta.size() < 9 ||
      u.size() < static_cast<std::size_t>(model.data.n_years)) {
    return c;
  }

  const int n_years = model.data.n_years;
  const int n_ages = model.data.n_ages;

  const double log_R0 = theta[0];
  const double log_M = theta[1];
  const double log_q = theta[2];
  const double log_Fbar = theta[3];
  const double sel50_raw = theta[4];
  const double log_sel_slope = theta[5];
  const double log_sigma_index = theta[6];
  const double log_sigma_catch = theta[7];
  const double log_sigma_rec = theta[8];

  const double R0 = std::exp(log_R0);
  const double M = std::exp(log_M);
  const double q = std::exp(log_q);
  const double Fbar = std::exp(log_Fbar);
  const double sel_slope = std::exp(log_sel_slope);

  const double sigma_index = 0.03 + std::exp(log_sigma_index);
  const double sigma_catch = 0.08 + std::exp(log_sigma_catch);
  const double sigma_rec = 0.05 + std::exp(log_sigma_rec);

  c.fixed_priors += normal_prior_nll_double(log_R0, std::log(900.0), 0.55);
  c.fixed_priors += normal_prior_nll_double(log_M, std::log(0.25), 0.22);
  c.fixed_priors += normal_prior_nll_double(log_q, std::log(0.15), 0.55);
  c.fixed_priors += normal_prior_nll_double(log_Fbar, std::log(0.18), 0.40);
  c.fixed_priors += normal_prior_nll_double(sel50_raw, 0.0, 1.50);
  c.fixed_priors +=
      normal_prior_nll_double(log_sel_slope, std::log(1.25), 0.75);
  c.fixed_priors +=
      normal_prior_nll_double(log_sigma_index, std::log(0.20), 0.75);
  c.fixed_priors +=
      normal_prior_nll_double(log_sigma_catch, std::log(0.18), 0.35);
  c.fixed_priors +=
      normal_prior_nll_double(log_sigma_rec, std::log(0.35), 0.75);

  double mean_rec_dev = 0.0;
  for (int y = 0; y < n_years; ++y) {
    mean_rec_dev += u[static_cast<std::size_t>(y)];
  }
  mean_rec_dev /= static_cast<double>(n_years);
  c.recruitment_mean_penalty += 0.5 * square_double(mean_rec_dev / 0.10);

  std::vector<double> numbers(static_cast<std::size_t>(n_ages),
                              R0 / static_cast<double>(n_ages));

  std::vector<double> selectivity(static_cast<std::size_t>(n_ages), 0.0);
  for (int a = 0; a < n_ages; ++a) {
    const double age = static_cast<double>(a + 1);
    const double sel50 = 3.0 + sel50_raw;
    selectivity[static_cast<std::size_t>(a)] =
        1.0 / (1.0 + std::exp(-(age - sel50) / sel_slope));
  }

  for (int y = 0; y < n_years; ++y) {
    const double rec_dev = u[static_cast<std::size_t>(y)];

    c.recruitment_prior +=
        0.5 * square_double(rec_dev / sigma_rec) + std::log(sigma_rec);

    c.recruitment_soft_stabilizer += 0.5 * square_double(rec_dev / 4.0);

    const double recruitment = R0 * std::exp(rec_dev);
    numbers[0] = recruitment;

    double vulnerable_biomass = 0.0;
    double total_catch = 0.0;

    for (int a = 0; a < n_ages; ++a) {
      const double age = static_cast<double>(a + 1);
      const double weight = (age * age * age) / 100.0;
      const double Na = numbers[static_cast<std::size_t>(a)];
      const double Sa = selectivity[static_cast<std::size_t>(a)];
      const double F_at_age = Fbar * Sa;
      const double Z = M + F_at_age;

      vulnerable_biomass += Na * weight * Sa;

      if (Z > 1.0e-12) {
        total_catch += Na * weight * (F_at_age / Z) * (1.0 - std::exp(-Z));
      }
    }

    const double pred_index = q * vulnerable_biomass + 1.0e-9;
    const double pred_catch = total_catch + 1.0e-9;

    const double obs_index = model.data.index_obs[static_cast<std::size_t>(y)];
    const double obs_catch = model.data.catch_obs[static_cast<std::size_t>(y)];

    c.index_likelihood += 0.5 * square_double((safe_log_double(obs_index) -
                                               safe_log_double(pred_index)) /
                                              sigma_index) +
                          std::log(sigma_index);

    c.catch_likelihood += 0.5 * square_double((safe_log_double(obs_catch) -
                                               safe_log_double(pred_catch)) /
                                              sigma_catch) +
                          std::log(sigma_catch);

    double vulnerable_number_sum = 0.0;
    for (int aa = 0; aa < n_ages; ++aa) {
      vulnerable_number_sum += numbers[static_cast<std::size_t>(aa)] *
                               selectivity[static_cast<std::size_t>(aa)];
    }

    for (int a = n_ages - 1; a >= 1; --a) {
      const double Sa_prev = selectivity[static_cast<std::size_t>(a - 1)];
      const double Z_prev = M + Fbar * Sa_prev;
      numbers[static_cast<std::size_t>(a)] =
          numbers[static_cast<std::size_t>(a - 1)] * std::exp(-Z_prev);
    }
  }

  return c;
}

// quadra_standalone_composition_proxy_v1
inline double
composition_proxy_from_age_diagnostic_basis(const CatchAtAgeLaplaceModel &model,
                                            const std::vector<double> &theta,
                                            const std::vector<double> &u) {
  if (theta.size() < 9 ||
      u.size() < static_cast<std::size_t>(model.data.n_years)) {
    return 0.0;
  }

  const int n_years = model.data.n_years;
  const int n_ages = model.data.n_ages;

  const double R0 = std::exp(theta[0]);
  const double M = std::exp(theta[1]);
  const double Fbar = std::exp(theta[3]);
  const double sel50_raw = theta[4];
  const double sel_slope = std::exp(theta[5]);

  std::vector<double> selectivity(static_cast<std::size_t>(n_ages), 0.0);
  for (int a = 0; a < n_ages; ++a) {
    const double age = static_cast<double>(a + 1);
    const double sel50 = 3.0 + sel50_raw;
    selectivity[static_cast<std::size_t>(a)] =
        1.0 / (1.0 + std::exp(-(age - sel50) / sel_slope));
  }

  std::vector<double> numbers(static_cast<std::size_t>(n_ages),
                              R0 / static_cast<double>(n_ages));

  double out = 0.0;

  for (int y = 0; y < n_years; ++y) {
    numbers[0] = R0 * std::exp(u[static_cast<std::size_t>(y)]);

    double vulnerable_number_sum = 0.0;
    for (int a = 0; a < n_ages; ++a) {
      vulnerable_number_sum += numbers[static_cast<std::size_t>(a)] *
                               selectivity[static_cast<std::size_t>(a)];
    }

    for (int a = 0; a < n_ages; ++a) {
      const double pred = (numbers[static_cast<std::size_t>(a)] *
                           selectivity[static_cast<std::size_t>(a)]) /
                          (vulnerable_number_sum + 1.0e-12);

      const double obs = model.data.age_comp_obs[static_cast<std::size_t>(y)]
                                                [static_cast<std::size_t>(a)];

      out += 0.5 * square_double((obs - pred) / 0.10);
    }

    for (int a = n_ages - 1; a >= 1; --a) {
      const double Sa_prev = selectivity[static_cast<std::size_t>(a - 1)];
      const double Z_prev = M + Fbar * Sa_prev;
      numbers[static_cast<std::size_t>(a)] =
          numbers[static_cast<std::size_t>(a - 1)] * std::exp(-Z_prev);
    }
  }

  return out;
}

template <typename Result>
void print_objective_decomposition_report(
    const CatchAtAgeLaplaceModel &model, const std::vector<double> &theta,
    const Result &result, const std::vector<double> &final_random_effects) {
  std::cout << std::endl;
  std::cout << "Objective decomposition report" << std::endl;
  {
    ObjectiveComponents c =
        evaluate_objective_components(model, theta, final_random_effects);

    c.composition_proxy = composition_proxy_from_age_diagnostic_basis(
        model, theta, final_random_effects);

    std::cout << std::setw(34) << "component" << std::setw(18) << "nll"
              << std::endl;

    std::cout << std::setw(34) << "fixed_priors" << std::setw(18)
              << c.fixed_priors << std::endl;
    std::cout << std::setw(34) << "recruitment_prior" << std::setw(18)
              << c.recruitment_prior << std::endl;
    std::cout << std::setw(34) << "recruitment_mean_penalty" << std::setw(18)
              << c.recruitment_mean_penalty << std::endl;
    std::cout << std::setw(34) << "recruitment_soft_stabilizer" << std::setw(18)
              << c.recruitment_soft_stabilizer << std::endl;
    std::cout << std::setw(34) << "index_likelihood" << std::setw(18)
              << c.index_likelihood << std::endl;
    std::cout << std::setw(34) << "catch_likelihood" << std::setw(18)
              << c.catch_likelihood << std::endl;
    std::cout << std::setw(34) << "composition_proxy" << std::setw(18)
              << c.composition_proxy << std::endl;
    std::cout << std::setw(34) << "component_total_joint" << std::setw(18)
              << c.total() << std::endl;
    std::cout << std::setw(34) << "reported_joint_objective" << std::setw(18)
              << result.joint_objective_m << std::endl;
    std::cout << std::setw(34) << "component_minus_reported" << std::setw(18)
              << c.total() - result.joint_objective_m << std::endl;

    std::cout
        << "Note: this is a diagnostic mirror, not the authoritative objective."
        << std::endl;
    std::cout << "      The authoritative check is model(fixed + random) "
                 "versus result.joint_objective_m."
              << std::endl;
    std::cout << "      Differences here mean the mirror uses different "
                 "bookkeeping than operator()."
              << std::endl;
  }
}

// quadra_big_example_age_comp_diagnostics_v1
template <typename Result>
void print_age_composition_diagnostics(
    const CatchAtAgeLaplaceModel &model, const std::vector<double> &theta,
    const Result &result, const std::vector<double> &final_random_effects) {
  std::cout << std::endl;
  std::cout << "Age-composition observed-vs-predicted diagnostics" << std::endl;

  if (theta.size() < 9) {
    std::cout
        << "age-composition diagnostics unavailable: expected 9 fixed effects"
        << std::endl;
    return;
  }
  {
    const auto &u = final_random_effects;

    const int n_years = model.data.n_years;
    const int n_ages = model.data.n_ages;

    if (u.size() < static_cast<std::size_t>(n_years)) {
      std::cout << "age-composition diagnostics unavailable: "
                   "final_random_effects shorter than n_years"
                << std::endl;
      return;
    }

    const double R0 = std::exp(theta[0]);
    const double M = std::exp(theta[1]);
    const double Fbar = std::exp(theta[3]);
    const double sel50_raw = theta[4];
    const double sel_slope = std::exp(theta[5]);

    std::vector<double> selectivity(static_cast<std::size_t>(n_ages), 0.0);
    for (int a = 0; a < n_ages; ++a) {
      const double age = static_cast<double>(a + 1);
      const double sel50 = 3.0 + sel50_raw;
      selectivity[static_cast<std::size_t>(a)] =
          1.0 / (1.0 + std::exp(-(age - sel50) / sel_slope));
    }

    std::vector<double> numbers(static_cast<std::size_t>(n_ages),
                                R0 / static_cast<double>(n_ages));

    double total_sse = 0.0;
    double total_abs = 0.0;
    double max_abs = 0.0;
    int n_resid = 0;

    std::cout << std::setw(6) << "year" << std::setw(12) << "row_sse"
              << std::setw(12) << "row_mae" << std::setw(12) << "row_max"
              << std::setw(14) << "obs_sum" << std::setw(14) << "pred_sum"
              << std::endl;

    for (int y = 0; y < n_years; ++y) {
      const double rec = R0 * std::exp(u[static_cast<std::size_t>(y)]);
      numbers[0] = rec;

      std::vector<double> pred(static_cast<std::size_t>(n_ages), 0.0);
      std::vector<double> obs(static_cast<std::size_t>(n_ages), 0.0);

      double pred_sum_raw = 0.0;
      double obs_sum = 0.0;

      // Use vulnerable numbers-at-age as the predicted composition basis.
      // This mirrors a common fishery age-composition observation model
      // more closely than dividing numbers-at-age by biomass.
      for (int a = 0; a < n_ages; ++a) {
        const double vulnerable_number =
            numbers[static_cast<std::size_t>(a)] *
            selectivity[static_cast<std::size_t>(a)];

        pred[static_cast<std::size_t>(a)] = vulnerable_number;
        pred_sum_raw += vulnerable_number;

        obs[static_cast<std::size_t>(a)] =
            model.data.age_comp_obs[static_cast<std::size_t>(y)]
                                   [static_cast<std::size_t>(a)];
        obs_sum += obs[static_cast<std::size_t>(a)];
      }

      for (int a = 0; a < n_ages; ++a) {
        pred[static_cast<std::size_t>(a)] =
            pred[static_cast<std::size_t>(a)] / (pred_sum_raw + 1.0e-12);
      }

      double pred_sum = 0.0;
      double row_sse = 0.0;
      double row_abs = 0.0;
      double row_max = 0.0;

      for (int a = 0; a < n_ages; ++a) {
        pred_sum += pred[static_cast<std::size_t>(a)];

        const double resid = obs[static_cast<std::size_t>(a)] -
                             pred[static_cast<std::size_t>(a)];

        const double abs_resid = std::abs(resid);
        row_sse += resid * resid;
        row_abs += abs_resid;
        row_max = std::max(row_max, abs_resid);

        total_sse += resid * resid;
        total_abs += abs_resid;
        max_abs = std::max(max_abs, abs_resid);
        ++n_resid;
      }

      std::cout << std::setw(6) << (y + 1) << std::setw(12) << row_sse
                << std::setw(12) << row_abs / static_cast<double>(n_ages)
                << std::setw(12) << row_max << std::setw(14) << obs_sum
                << std::setw(14) << pred_sum << std::endl;

      if (y < 5 || y == n_years - 1) {
        std::cout << "  age detail:";
        for (int a = 0; a < n_ages; ++a) {
          const double resid = obs[static_cast<std::size_t>(a)] -
                               pred[static_cast<std::size_t>(a)];

          std::cout << " a" << (a + 1)
                    << "(o=" << obs[static_cast<std::size_t>(a)]
                    << ",p=" << pred[static_cast<std::size_t>(a)]
                    << ",r=" << resid << ")";
        }
        std::cout << std::endl;
      }

      for (int a = n_ages - 1; a >= 1; --a) {
        const double Sa_prev = selectivity[static_cast<std::size_t>(a - 1)];
        const double Z_prev = M + Fbar * Sa_prev;
        numbers[static_cast<std::size_t>(a)] =
            numbers[static_cast<std::size_t>(a - 1)] * std::exp(-Z_prev);
      }
    }

    std::cout << std::endl;
    std::cout << "Age-composition residual summary" << std::endl;
    std::cout << "total_sse: " << total_sse << std::endl;
    std::cout << "mean_abs_residual: "
              << (n_resid > 0 ? total_abs / static_cast<double>(n_resid) : 0.0)
              << std::endl;
    std::cout << "max_abs_residual: " << max_abs << std::endl;
    std::cout << "proxy_nll_at_sd_0.10: " << 0.5 * total_sse / (0.10 * 0.10)
              << std::endl;
    std::cout
        << "Note: proxy_nll_at_sd_0.10 should be compared to composition_proxy."
        << std::endl;
    std::cout << "      Large discrepancies indicate the diagnostic mirror and "
                 "model objective use different predictions."
              << std::endl;
  }
}

// quadra_index_catch_diagnostics_v1
template <typename Result>
void print_index_catch_diagnostics(
    const CatchAtAgeLaplaceModel &model, const std::vector<double> &theta,
    const Result &result, const std::vector<double> &final_random_effects) {
  std::cout << std::endl;
  std::cout << "Index and catch observed-vs-predicted diagnostics" << std::endl;

  if (theta.size() < 9) {
    std::cout << "index/catch diagnostics unavailable: expected 9 fixed effects"
              << std::endl;
    return;
  }
  {
    const auto &u = final_random_effects;

    const int n_years = model.data.n_years;
    const int n_ages = model.data.n_ages;

    if (u.size() < static_cast<std::size_t>(n_years)) {
      std::cout << "index/catch diagnostics unavailable: final_random_effects "
                   "shorter than n_years"
                << std::endl;
      return;
    }

    const double R0 = std::exp(theta[0]);
    const double M = std::exp(theta[1]);
    const double q = std::exp(theta[2]);
    const double Fbar = std::exp(theta[3]);
    const double sel50_raw = theta[4];
    const double sel_slope = std::exp(theta[5]);
    const double log_sigma_index = theta[6];
    const double log_sigma_catch = theta[7];

    const double sigma_index = 0.03 + std::exp(log_sigma_index);
    const double sigma_catch = 0.08 + std::exp(log_sigma_catch);

    std::vector<double> selectivity(static_cast<std::size_t>(n_ages), 0.0);
    for (int a = 0; a < n_ages; ++a) {
      const double age = static_cast<double>(a + 1);
      const double sel50 = 3.0 + sel50_raw;
      selectivity[static_cast<std::size_t>(a)] =
          1.0 / (1.0 + std::exp(-(age - sel50) / sel_slope));
    }

    std::vector<double> numbers(static_cast<std::size_t>(n_ages),
                                R0 / static_cast<double>(n_ages));

    double total_index_nll = 0.0;
    double total_catch_nll = 0.0;
    double total_index_sq = 0.0;
    double total_catch_sq = 0.0;
    double max_abs_index_log_resid = 0.0;
    double max_abs_catch_log_resid = 0.0;

    std::cout << "sigma_index_used: " << sigma_index << std::endl;
    std::cout << "sigma_catch_used: " << sigma_catch << std::endl;

    std::cout << std::endl;
    std::cout << std::setw(6) << "year" << std::setw(14) << "obs_index"
              << std::setw(14) << "pred_index" << std::setw(14) << "idx_logres"
              << std::setw(14) << "idx_nll" << std::setw(14) << "obs_catch"
              << std::setw(14) << "pred_catch" << std::setw(14) << "cat_logres"
              << std::setw(14) << "cat_nll" << std::endl;

    for (int y = 0; y < n_years; ++y) {
      numbers[0] = R0 * std::exp(u[static_cast<std::size_t>(y)]);

      double vulnerable_biomass = 0.0;
      double total_catch = 0.0;

      for (int a = 0; a < n_ages; ++a) {
        const double age = static_cast<double>(a + 1);
        const double weight = (age * age * age) / 100.0;
        const double Na = numbers[static_cast<std::size_t>(a)];
        const double Sa = selectivity[static_cast<std::size_t>(a)];
        const double F_at_age = Fbar * Sa;
        const double Z = M + F_at_age;

        vulnerable_biomass += Na * weight * Sa;

        if (Z > 1.0e-12) {
          total_catch += Na * weight * (F_at_age / Z) * (1.0 - std::exp(-Z));
        }
      }

      const double pred_index = q * vulnerable_biomass + 1.0e-9;
      const double pred_catch = total_catch + 1.0e-9;

      const double obs_index =
          model.data.index_obs[static_cast<std::size_t>(y)];
      const double obs_catch =
          model.data.catch_obs[static_cast<std::size_t>(y)];

      const double index_log_resid =
          safe_log_double(obs_index) - safe_log_double(pred_index);
      const double catch_log_resid =
          safe_log_double(obs_catch) - safe_log_double(pred_catch);

      const double index_nll =
          0.5 * square_double(index_log_resid / sigma_index) +
          std::log(sigma_index);
      const double catch_nll =
          0.5 * square_double(catch_log_resid / sigma_catch) +
          std::log(sigma_catch);

      total_index_nll += index_nll;
      total_catch_nll += catch_nll;
      total_index_sq += index_log_resid * index_log_resid;
      total_catch_sq += catch_log_resid * catch_log_resid;
      max_abs_index_log_resid =
          std::max(max_abs_index_log_resid, std::abs(index_log_resid));
      max_abs_catch_log_resid =
          std::max(max_abs_catch_log_resid, std::abs(catch_log_resid));

      std::cout << std::setw(6) << (y + 1) << std::setw(14) << obs_index
                << std::setw(14) << pred_index << std::setw(14)
                << index_log_resid << std::setw(14) << index_nll
                << std::setw(14) << obs_catch << std::setw(14) << pred_catch
                << std::setw(14) << catch_log_resid << std::setw(14)
                << catch_nll << std::endl;

      for (int a = n_ages - 1; a >= 1; --a) {
        const double Sa_prev = selectivity[static_cast<std::size_t>(a - 1)];
        const double Z_prev = M + Fbar * Sa_prev;
        numbers[static_cast<std::size_t>(a)] =
            numbers[static_cast<std::size_t>(a - 1)] * std::exp(-Z_prev);
      }
    }

    std::cout << std::endl;
    std::cout << "Index/catch residual summary" << std::endl;
    std::cout << "index_nll_total: " << total_index_nll << std::endl;
    std::cout << "catch_nll_total: " << total_catch_nll << std::endl;
    std::cout << "index_logres_rmse: "
              << std::sqrt(total_index_sq / static_cast<double>(n_years))
              << std::endl;
    std::cout << "catch_logres_rmse: "
              << std::sqrt(total_catch_sq / static_cast<double>(n_years))
              << std::endl;
    std::cout << "index_max_abs_logres: " << max_abs_index_log_resid
              << std::endl;
    std::cout << "catch_max_abs_logres: " << max_abs_catch_log_resid
              << std::endl;
    std::cout << "Note: index_nll_total and catch_nll_total should match the "
                 "decomposition report."
              << std::endl;
  }
}

} // namespace example

// quadra_actual_objective_path_decomposition_v1
struct ActualObjectivePathComponents {
  double fixed_priors = 0.0;
  double recruitment_prior = 0.0;
  double recruitment_mean_penalty = 0.0;
  double recruitment_soft_stabilizer = 0.0;
  double index_likelihood = 0.0;
  double catch_likelihood = 0.0;
  double composition_likelihood = 0.0;

  double total() const {
    return fixed_priors + recruitment_prior + recruitment_mean_penalty +
           recruitment_soft_stabilizer + index_likelihood + catch_likelihood +
           composition_likelihood;
  }
};

inline ActualObjectivePathComponents evaluate_actual_objective_path_components(
    const example::CatchAtAgeLaplaceModel &model,
    const std::vector<double> &theta, const std::vector<double> &u) {
  ActualObjectivePathComponents c;

  if (theta.size() < 9 ||
      u.size() < static_cast<std::size_t>(model.data.n_years)) {
    return c;
  }

  const int n_years = model.data.n_years;
  const int n_ages = model.data.n_ages;

  const double log_R0 = theta[0];
  const double log_M = theta[1];
  const double log_q = theta[2];
  const double log_Fbar = theta[3];
  const double sel50_raw = theta[4];
  const double log_sel_slope = theta[5];
  const double log_sigma_index = theta[6];
  const double log_sigma_catch = theta[7];
  const double log_sigma_rec = theta[8];

  const double R0 = std::exp(log_R0);
  const double M = std::exp(log_M);
  const double q = std::exp(log_q);
  const double Fbar = std::exp(log_Fbar);
  const double sel_slope = std::exp(log_sel_slope);

  const double sigma_index = 0.03 + std::exp(log_sigma_index);
  const double sigma_catch = 0.08 + std::exp(log_sigma_catch);
  const double sigma_rec = 0.05 + std::exp(log_sigma_rec);

  c.fixed_priors +=
      example::normal_prior_nll_double(log_R0, std::log(900.0), 0.55);
  c.fixed_priors +=
      example::normal_prior_nll_double(log_M, std::log(0.25), 0.22);
  c.fixed_priors +=
      example::normal_prior_nll_double(log_q, std::log(0.15), 0.55);
  c.fixed_priors +=
      example::normal_prior_nll_double(log_Fbar, std::log(0.18), 0.40);
  c.fixed_priors += example::normal_prior_nll_double(sel50_raw, 0.0, 1.50);
  c.fixed_priors +=
      example::normal_prior_nll_double(log_sel_slope, std::log(1.25), 0.75);
  c.fixed_priors +=
      example::normal_prior_nll_double(log_sigma_index, std::log(0.20), 0.75);
  c.fixed_priors +=
      example::normal_prior_nll_double(log_sigma_catch, std::log(0.18), 0.35);
  c.fixed_priors +=
      example::normal_prior_nll_double(log_sigma_rec, std::log(0.35), 0.75);

  double mean_rec_dev = 0.0;
  for (int y = 0; y < n_years; ++y) {
    mean_rec_dev += u[static_cast<std::size_t>(y)];
  }
  mean_rec_dev /= static_cast<double>(n_years);
  c.recruitment_mean_penalty +=
      0.5 * example::square_double(mean_rec_dev / 0.10);

  std::vector<double> numbers(static_cast<std::size_t>(n_ages),
                              R0 / static_cast<double>(n_ages));

  std::vector<double> selectivity(static_cast<std::size_t>(n_ages), 0.0);
  for (int a = 0; a < n_ages; ++a) {
    const double age = static_cast<double>(a + 1);
    const double sel50 = 3.0 + sel50_raw;
    selectivity[static_cast<std::size_t>(a)] =
        1.0 / (1.0 + std::exp(-(age - sel50) / sel_slope));
  }

  for (int y = 0; y < n_years; ++y) {
    const double rec_dev = u[static_cast<std::size_t>(y)];

    c.recruitment_prior +=
        0.5 * example::square_double(rec_dev / sigma_rec) + std::log(sigma_rec);

    c.recruitment_soft_stabilizer +=
        0.5 * example::square_double(rec_dev / 4.0);

    numbers[0] = R0 * std::exp(rec_dev);

    double vulnerable_biomass = 0.0;
    double total_catch = 0.0;

    for (int a = 0; a < n_ages; ++a) {
      const double age = static_cast<double>(a + 1);
      const double weight = (age * age * age) / 100.0;
      const double Na = numbers[static_cast<std::size_t>(a)];
      const double Sa = selectivity[static_cast<std::size_t>(a)];
      const double F_at_age = Fbar * Sa;
      const double Z = M + F_at_age;

      vulnerable_biomass += Na * weight * Sa;

      if (Z > 1.0e-12) {
        total_catch += Na * weight * (F_at_age / Z) * (1.0 - std::exp(-Z));
      }
    }

    const double pred_index = q * vulnerable_biomass + 1.0e-9;
    const double pred_catch = total_catch + 1.0e-9;

    const double obs_index = model.data.index_obs[static_cast<std::size_t>(y)];
    const double obs_catch = model.data.catch_obs[static_cast<std::size_t>(y)];

    c.index_likelihood +=
        0.5 * example::square_double((example::safe_log_double(obs_index) -
                                      example::safe_log_double(pred_index)) /
                                     sigma_index) +
        std::log(sigma_index);

    c.catch_likelihood +=
        0.5 * example::square_double((example::safe_log_double(obs_catch) -
                                      example::safe_log_double(pred_catch)) /
                                     sigma_catch) +
        std::log(sigma_catch);

    double vulnerable_number_sum = 0.0;
    for (int a = 0; a < n_ages; ++a) {
      vulnerable_number_sum += numbers[static_cast<std::size_t>(a)] *
                               selectivity[static_cast<std::size_t>(a)];
    }

    for (int a = 0; a < n_ages; ++a) {
      const double pred = (numbers[static_cast<std::size_t>(a)] *
                           selectivity[static_cast<std::size_t>(a)]) /
                          (vulnerable_number_sum + 1.0e-12);

      const double obs = model.data.age_comp_obs[static_cast<std::size_t>(y)]
                                                [static_cast<std::size_t>(a)];

      c.composition_likelihood +=
          0.5 * example::square_double((obs - pred) / 0.10);
    }

    for (int a = n_ages - 1; a >= 1; --a) {
      const double Sa_prev = selectivity[static_cast<std::size_t>(a - 1)];
      const double Z_prev = M + Fbar * Sa_prev;
      numbers[static_cast<std::size_t>(a)] =
          numbers[static_cast<std::size_t>(a - 1)] * std::exp(-Z_prev);
    }
  }

  return c;
}

template <typename Result>
void print_actual_objective_path_decomposition(
    const example::CatchAtAgeLaplaceModel &model,
    const std::vector<double> &theta, const Result &result,
    const std::vector<double> &final_random_effects) {
  std::cout << std::endl;
  std::cout << "Diagnostic objective-path decomposition" << std::endl;
  {
    const auto c = evaluate_actual_objective_path_components(
        model, theta, final_random_effects);

    std::cout << std::setw(34) << "component" << std::setw(18) << "nll"
              << std::endl;

    std::cout << std::setw(34) << "fixed_priors" << std::setw(18)
              << c.fixed_priors << std::endl;
    std::cout << std::setw(34) << "recruitment_prior" << std::setw(18)
              << c.recruitment_prior << std::endl;
    std::cout << std::setw(34) << "recruitment_mean_penalty" << std::setw(18)
              << c.recruitment_mean_penalty << std::endl;
    std::cout << std::setw(34) << "recruitment_soft_stabilizer" << std::setw(18)
              << c.recruitment_soft_stabilizer << std::endl;
    std::cout << std::setw(34) << "index_likelihood" << std::setw(18)
              << c.index_likelihood << std::endl;
    std::cout << std::setw(34) << "catch_likelihood" << std::setw(18)
              << c.catch_likelihood << std::endl;
    std::cout << std::setw(34) << "composition_likelihood" << std::setw(18)
              << c.composition_likelihood << std::endl;
    std::cout << std::setw(34) << "diagnostic_path_total" << std::setw(18)
              << c.total() << std::endl;
    std::cout << std::setw(34) << "reported_joint_objective" << std::setw(18)
              << result.joint_objective_m << std::endl;
    std::cout << std::setw(34) << "diagnostic_path_minus_reported"
              << std::setw(18) << c.total() - result.joint_objective_m
              << std::endl;
  }
}

namespace example {

inline quadra::ParameterVector make_big_laplace_parameter_vector() {
  quadra::ParameterVector param_vector;

  param_vector.add(quadra::Parameter(
      "log_R0", std::log(900.0), quadra::ParameterTransform::Identity, false));
  param_vector.add(quadra::Parameter(
      "log_M", std::log(0.25), quadra::ParameterTransform::Identity, false));
  param_vector.add(quadra::Parameter(
      "log_q", std::log(0.15), quadra::ParameterTransform::Identity, false));
  param_vector.add(quadra::Parameter(
      "log_Fbar", std::log(0.18), quadra::ParameterTransform::Identity, false));
  param_vector.add(quadra::Parameter(
      "sel50_raw", 0.0, quadra::ParameterTransform::Identity, false));
  param_vector.add(quadra::Parameter("log_sel_slope", std::log(1.25),
                                     quadra::ParameterTransform::Identity,
                                     false));
  param_vector.add(quadra::Parameter("log_sigma_index", std::log(0.20),
                                     quadra::ParameterTransform::Identity,
                                     false));
  param_vector.add(quadra::Parameter("log_sigma_catch", std::log(0.18),
                                     quadra::ParameterTransform::Identity,
                                     false));
  param_vector.add(quadra::Parameter("log_sigma_rec", std::log(0.35),
                                     quadra::ParameterTransform::Identity,
                                     false));

  for (int y = 0; y < 30; ++y) {
    param_vector.add(quadra::Parameter("rec_dev_" + std::to_string(y + 1), 0.0,
                                       quadra::ParameterTransform::Identity,
                                       true));
  }

  return param_vector;
}

} // namespace example
