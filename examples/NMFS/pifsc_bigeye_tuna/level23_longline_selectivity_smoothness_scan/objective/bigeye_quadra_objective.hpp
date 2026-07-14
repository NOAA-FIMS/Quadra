#pragma once

#include "../quadra/bigeye_age_structured.hpp"

#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace pifsc_bigeye_tuna {

struct FleetObservation {
  int year = 0;
  std::string fleet;
  double catch_mt = 0.0;
  double index = 0.0;
  std::array<double, kAges> age_comp{};
};

template <class T> T exp_t(const T &x) {
  using std::exp;
  return exp(x);
}
template <class T> T log_t(const T &x) {
  using std::log;
  return log(x);
}
template <class T> T invlogit_t(const T &x) {
  return T(1.0) / (T(1.0) + exp_t(-x));
}
template <class T> T max_t(const T &x, double floor) {
  return x > T(floor) ? x : T(floor);
}
template <class T> T square_t(const T &x) { return x * x; }

template <class T>
T logistic_selectivity_t(const T &age, const T &a50, const T &slope) {
  return T(1.0) / (T(1.0) + exp_t(-slope * (age - a50)));
}

template <class T>
T age_comp_nll(const std::array<double, kAges> &observed,
               const std::array<T, kAges> &predicted, double effective_n,
               double floor = 1.0e-12) {
  T nll = T(0.0);
  for (int a = 0; a < kAges; ++a) {
    const auto i = static_cast<std::size_t>(a);
    const double obs = std::max(observed[i], 0.0);
    if (obs > 0.0)
      nll = nll - T(effective_n * obs) * log_t(max_t(predicted[i], floor));
  }
  return nll;
}

template <class T> std::array<T, kAges> fixed_purse_seine_age_selectivity() {
  const std::array<double, kAges> raw = {1.00, 1.00, 0.85, 0.55, 0.25,
                                         0.10, 0.04, 0.02, 0.01, 0.005};
  std::array<T, kAges> out{};
  for (int a = 0; a < kAges; ++a)
    out[static_cast<std::size_t>(a)] = T(raw[static_cast<std::size_t>(a)]);
  return out;
}

class BigeyeQuadraObjective {
public:
  explicit BigeyeQuadraObjective(std::vector<FleetObservation> observations)
      : observations_(std::move(observations)) {
    if (observations_.empty())
      throw std::runtime_error("BigeyeQuadraObjective requires observations");
  }

  template <class T> T operator()(const std::vector<T> &par) const {
    constexpr int kBaseFixed = 3;
    constexpr int kLonglineSelOffset = kBaseFixed;
    constexpr int kLonglineSelDevs = kAges;
    constexpr int kInitialDevs = kAges;
    constexpr int kInitialDevOffset = kLonglineSelOffset + kLonglineSelDevs;
    constexpr int kPurseSeineSelOffset = kInitialDevOffset + kInitialDevs;
    constexpr int kPurseSeineSelDevs = kAges;
    constexpr int kRecruitmentOffset =
        kPurseSeineSelOffset + kPurseSeineSelDevs;

    if (par.size() < static_cast<std::size_t>(kRecruitmentOffset) + n_years())
      throw std::runtime_error(
          "Level 23 expected 3 base fixed effects, longline age selectivity "
          "logits, initial number deviations, purse-seine age selectivity "
          "logits, plus recruitment deviations");

    const T log_r0 = par[0];
    const T log_fbar = par[1];
    const T log_q_purse_seine = par[2];

    const T log_m = T(std::log(0.45));             // fixed M=0.45 anchor
    const T log_q_longline = T(std::log(0.00005)); // fixed q anchor

    const T r0 = exp_t(log_r0);
    const T m = exp_t(log_m);
    const T fbar = exp_t(log_fbar);
    const T q_longline = exp_t(log_q_longline);
    const T q_purse_seine = exp_t(log_q_purse_seine);

    const T sigma_log_index = T(0.20);
    const T sigma_log_catch = T(0.15);
    const T sigma_rec_dev = T(0.35);
    const T steepness = T(0.75);
    const T sigma_init_dev = T(0.75);
    const double age_comp_effective_n = 30.0;
    const double min_positive = 1.0e-12;

    const auto weight = default_weight_at_age();
    const auto maturity = default_maturity_at_age();

    const std::array<double, kAges> ll_template = {
        0.001, 0.005, 0.03, 0.12, 0.35, 0.70, 0.95, 0.90, 0.60, 0.12};

    std::array<T, kAges> sel_longline{};
    for (int a = 0; a < kAges; ++a) {
      sel_longline[static_cast<std::size_t>(a)] =
          invlogit_t(par[kLonglineSelOffset + a]);
    }

    const std::array<double, kAges> ps_template = {
        0.20, 0.90, 1.00, 0.80, 0.45, 0.20, 0.08, 0.03, 0.015, 0.005};

    std::array<T, kAges> sel_purse_seine{};
    for (int a = 0; a < kAges; ++a) {
      sel_purse_seine[static_cast<std::size_t>(a)] =
          invlogit_t(par[kPurseSeineSelOffset + a]);
    }

    T nll = T(0.0);
    auto normal_prior = [](const T &x, double mean, double sd) {
      const T z = (x - T(mean)) / T(sd);
      return T(0.5) * z * z;
    };

    nll = nll + normal_prior(log_r0, std::log(1200.0), 1.0);
    nll = nll + normal_prior(log_fbar, std::log(0.025), 0.75);
    nll = nll + normal_prior(log_q_purse_seine, std::log(0.00005), 1.0);

#ifdef BIGEYE_LL_SEL_SIGMA
    const T sigma_ll_sel_dev = T(BIGEYE_LL_SEL_SIGMA);
#else
    double sigma_ll_sel_dev_raw = 1.75;
    if (const char *raw_sigma_ll = std::getenv("LEVEL23_LL_SEL_SIGMA")) {
      try {
        sigma_ll_sel_dev_raw = std::stod(raw_sigma_ll);
      } catch (...) {
        sigma_ll_sel_dev_raw = 1.75;
      }
    }
    const T sigma_ll_sel_dev = T(sigma_ll_sel_dev_raw);
#endif
    for (int a = 0; a < kAges; ++a) {
      const auto i = static_cast<std::size_t>(a);
      const double p0 = std::min(0.999, std::max(0.001, ll_template[i]));
      const T prior_logit = T(std::log(p0 / (1.0 - p0)));
      const T raw_logit = par[kLonglineSelOffset + a];
      nll =
          nll + T(0.5) * square_t((raw_logit - prior_logit) / sigma_ll_sel_dev);
    }

    double ll_smooth_lambda_raw = 0.0;
    if (const char *raw_ll_smooth = std::getenv("LEVEL23_LL_SMOOTH_LAMBDA")) {
      try {
        ll_smooth_lambda_raw = std::stod(raw_ll_smooth);
      } catch (...) {
        ll_smooth_lambda_raw = 0.0;
      }
    }
    const T ll_smooth_lambda = T(ll_smooth_lambda_raw);
    for (int a = 1; a < kAges - 1; ++a) {
      const T left = par[kLonglineSelOffset + a - 1];
      const T center = par[kLonglineSelOffset + a];
      const T right = par[kLonglineSelOffset + a + 1];
      const T second_diff = left - T(2.0) * center + right;
      nll = nll + T(0.5) * ll_smooth_lambda * second_diff * second_diff;
    }

    const T sigma_ps_sel_dev = T(1.0);
    for (int a = 0; a < kAges; ++a) {
      const auto i = static_cast<std::size_t>(a);
      const double p0 = std::min(0.999, std::max(0.001, ps_template[i]));
      const T prior_logit = T(std::log(p0 / (1.0 - p0)));
      const T raw_logit = par[kPurseSeineSelOffset + a];
      nll =
          nll + T(0.5) * square_t((raw_logit - prior_logit) / sigma_ps_sel_dev);
    }

    std::array<T, kAges> n{};
    n[0] = r0;
    for (int a = 1; a < kAges; ++a)
      n[static_cast<std::size_t>(a)] =
          n[static_cast<std::size_t>(a - 1)] * exp_t(-m);
    n[static_cast<std::size_t>(kAges - 1)] =
        n[static_cast<std::size_t>(kAges - 1)] / (T(1.0) - exp_t(-m));

    T unfished_spawning_biomass = T(0.0);
    for (int a = 0; a < kAges; ++a) {
      const auto i = static_cast<std::size_t>(a);
      unfished_spawning_biomass =
          unfished_spawning_biomass + n[i] * T(weight[i]) * T(maturity[i]);
    }
    const T phi0 = unfished_spawning_biomass / max_t(r0, min_positive);

    for (int a = 0; a < kAges; ++a) {
      const T init_dev = par[kInitialDevOffset + a];
      nll = nll + T(0.5) * square_t(init_dev / sigma_init_dev);
      n[static_cast<std::size_t>(a)] =
          n[static_cast<std::size_t>(a)] * exp_t(init_dev);
    }

    const auto years = unique_years();

    for (std::size_t t = 0; t < years.size(); ++t) {
      const T rec_dev = par[kRecruitmentOffset + t];
      nll = nll + T(0.5) * square_t(rec_dev / sigma_rec_dev);

      std::array<T, kAges> catch_at_age_longline{};
      std::array<T, kAges> catch_at_age_purse_seine{};
      T longline_catch_hat = T(0.0);
      T purse_seine_catch_hat = T(0.0);

      for (int a = 0; a < kAges; ++a) {
        const auto i = static_cast<std::size_t>(a);
        const T f_longline_a = fbar * sel_longline[i];
        const T f_purse_seine_a = fbar * sel_purse_seine[i];
        const T f_total_a = f_longline_a + f_purse_seine_a;
        const T z_a = m + f_total_a;
        const T common_catch_factor =
            n[i] * (T(1.0) - exp_t(-z_a)) / max_t(z_a, min_positive);

        catch_at_age_longline[i] = common_catch_factor * f_longline_a;
        catch_at_age_purse_seine[i] = common_catch_factor * f_purse_seine_a;

        longline_catch_hat =
            longline_catch_hat + catch_at_age_longline[i] * T(weight[i]);
        purse_seine_catch_hat =
            purse_seine_catch_hat + catch_at_age_purse_seine[i] * T(weight[i]);
      }

      for (const auto &obs : observations_) {
        if (obs.year != years[t])
          continue;

        const bool is_longline = obs.fleet == "longline";
        const auto &sel = is_longline ? sel_longline : sel_purse_seine;
        const T fleet_q = is_longline ? q_longline : q_purse_seine;
        const auto &fleet_catch_at_age =
            is_longline ? catch_at_age_longline : catch_at_age_purse_seine;

        T vulnerable_biomass = T(0.0);
        T selected_numbers_sum = T(0.0);
        std::array<T, kAges> pred_age_comp{};

        for (int a = 0; a < kAges; ++a) {
          const auto i = static_cast<std::size_t>(a);
          vulnerable_biomass =
              vulnerable_biomass + n[i] * T(weight[i]) * sel[i];
          pred_age_comp[i] = fleet_catch_at_age[i];
          selected_numbers_sum = selected_numbers_sum + pred_age_comp[i];
        }

        const T index_hat = fleet_q * vulnerable_biomass;
        if (obs.index > 0.0) {
          const T z =
              (log_t(T(obs.index)) - log_t(max_t(index_hat, min_positive))) /
              sigma_log_index;
          nll = nll + T(0.5) * square_t(z);
        }

        const T catch_hat =
            is_longline ? longline_catch_hat : purse_seine_catch_hat;
        if (obs.catch_mt > 0.0) {
          const T z =
              (log_t(T(obs.catch_mt)) - log_t(max_t(catch_hat, min_positive))) /
              sigma_log_catch;
          nll = nll + T(0.5) * square_t(z);
        }

        for (int a = 0; a < kAges; ++a) {
          const auto i = static_cast<std::size_t>(a);
          pred_age_comp[i] =
              pred_age_comp[i] / max_t(selected_numbers_sum, min_positive);
        }

        nll = nll + age_comp_nll(obs.age_comp, pred_age_comp,
                                 age_comp_effective_n, min_positive);
      }

      T spawning_biomass = T(0.0);
      for (int a = 0; a < kAges; ++a) {
        const auto i = static_cast<std::size_t>(a);
        spawning_biomass =
            spawning_biomass + n[i] * T(weight[i]) * T(maturity[i]);
      }

      const T expected_recruitment =
          (T(0.8) * r0 * steepness * spawning_biomass) /
          max_t(T(0.2) * r0 * phi0 * (T(1.0) - steepness) +
                    spawning_biomass * (steepness - T(0.2)),
                min_positive);

      std::array<T, kAges> next{};
      next[0] = expected_recruitment * exp_t(rec_dev);

      for (int a = 1; a < kAges; ++a) {
        const auto prev = static_cast<std::size_t>(a - 1);
        const T total_sel_prev = sel_longline[prev] + sel_purse_seine[prev];
        const T z_prev = m + fbar * total_sel_prev;
        next[static_cast<std::size_t>(a)] = n[prev] * exp_t(-z_prev);
      }

      const auto last = static_cast<std::size_t>(kAges - 1);
      const T total_sel_last = sel_longline[last] + sel_purse_seine[last];
      const T z_last = m + fbar * total_sel_last;
      next[last] = next[last] + n[last] * exp_t(-z_last);
      n = next;
    }

    return nll;
  }

  std::size_t n_years() const { return unique_years().size(); }

  std::vector<int> unique_years() const {
    std::vector<int> years;
    for (const auto &obs : observations_) {
      if (years.empty() || years.back() != obs.year)
        years.push_back(obs.year);
    }
    return years;
  }

  const std::vector<FleetObservation> &fleet_observations() const {
    return observations_;
  }

private:
  std::vector<FleetObservation> observations_;
};

} // namespace pifsc_bigeye_tuna

using pifsc_bigeye_tuna::BigeyeQuadraObjective;
using pifsc_bigeye_tuna::FleetObservation;
