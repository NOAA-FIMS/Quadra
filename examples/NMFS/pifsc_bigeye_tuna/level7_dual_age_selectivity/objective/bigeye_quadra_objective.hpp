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

template <class T> T max_t(const T &x, double floor) {
  return x > T(floor) ? x : T(floor);
}

template <class T> T square_t(const T &x) { return x * x; }

template <class T>
T age_comp_nll(const std::array<double, kAges> &observed,
               const std::array<T, kAges> &predicted, double effective_n,
               double floor = 1.0e-12) {
  T nll = T(0.0);
  for (int a = 0; a < kAges; ++a) {
    const auto i = static_cast<std::size_t>(a);
    const double obs = std::max(observed[i], 0.0);
    if (obs > 0.0) {
      nll = nll - T(effective_n * obs) * log_t(max_t(predicted[i], floor));
    }
  }
  return nll;
}

template <class T> std::array<T, kAges> fixed_longline_age_selectivity() {
  // Older-fish longline vulnerability pattern.
  // Normalized to max = 1.0. This deliberately removes the Level 6
  // two-parameter longline logistic curve to test whether recruitment
  // persistence was compensating for remaining selectivity misspecification.
  const std::array<double, kAges> raw = {0.01, 0.03, 0.08, 0.18, 0.40,
                                         0.70, 0.95, 1.00, 0.95, 0.85};

  std::array<T, kAges> out{};
  for (int a = 0; a < kAges; ++a)
    out[static_cast<std::size_t>(a)] = T(raw[static_cast<std::size_t>(a)]);
  return out;
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
    if (observations_.empty()) {
      throw std::runtime_error("BigeyeQuadraObjective requires observations");
    }
  }

  template <class T> T operator()(const std::vector<T> &par) const {
    if (par.size() < 4 + n_years()) {
      throw std::runtime_error("BigeyeQuadraObjective expected 4 fixed effects "
                               "plus recruitment deviations");
    }

    const T log_r0 = par[0];
    const T log_fbar = par[1];
    const T log_q_longline = par[2];
    const T log_q_purse_seine = par[3];

    const T r0 = exp_t(log_r0);
    const T m = T(0.18);
    const T fbar = exp_t(log_fbar);
    const T q_longline = exp_t(log_q_longline);
    const T q_purse_seine = exp_t(log_q_purse_seine);

    const T sigma_log_index = T(0.20);
    const T sigma_log_catch = T(0.15);
    const T sigma_rec_dev = T(0.35);
    const double age_comp_effective_n = 30.0;
    const double min_positive = 1.0e-12;

    const auto weight = default_weight_at_age();
    const auto sel_longline = fixed_longline_age_selectivity<T>();
    const auto sel_purse_seine = fixed_purse_seine_age_selectivity<T>();

    std::array<T, kAges> n{};
    n[0] = r0;
    for (int a = 1; a < kAges; ++a) {
      n[static_cast<std::size_t>(a)] =
          n[static_cast<std::size_t>(a - 1)] * exp_t(-m);
    }
    n[static_cast<std::size_t>(kAges - 1)] =
        n[static_cast<std::size_t>(kAges - 1)] / (T(1.0) - exp_t(-m));

    T nll = T(0.0);

    auto normal_prior = [](const T &x, double mean, double sd) {
      const T z = (x - T(mean)) / T(sd);
      return T(0.5) * z * z;
    };

    nll = nll + normal_prior(log_r0, std::log(1200.0), 1.0);
    nll = nll + normal_prior(log_fbar, std::log(0.025), 0.75);
    nll = nll + normal_prior(log_q_longline, std::log(0.00005), 1.0);
    nll = nll + normal_prior(log_q_purse_seine, std::log(0.00005), 1.0);

    const auto years = unique_years();

    for (std::size_t t = 0; t < years.size(); ++t) {
      const T rec_dev = par[4 + t];
      nll = nll + T(0.5) * square_t(rec_dev / sigma_rec_dev);

      T total_catch_hat = T(0.0);

      for (int a = 0; a < kAges; ++a) {
        const auto i = static_cast<std::size_t>(a);
        const T avg_sel = T(0.5) * (sel_longline[i] + sel_purse_seine[i]);
        const T f_a = fbar * avg_sel;
        const T z_a = m + f_a;
        const T harvest_rate = (f_a / z_a) * (T(1.0) - exp_t(-z_a));
        total_catch_hat = total_catch_hat + n[i] * T(weight[i]) * harvest_rate;
      }

      for (const auto &obs : observations_) {
        if (obs.year != years[t]) {
          continue;
        }

        const bool is_longline = obs.fleet == "longline";
        const auto &sel = is_longline ? sel_longline : sel_purse_seine;
        const T fleet_q = is_longline ? q_longline : q_purse_seine;

        T vulnerable_biomass = T(0.0);
        T selected_numbers_sum = T(0.0);
        std::array<T, kAges> pred_age_comp{};

        for (int a = 0; a < kAges; ++a) {
          const auto i = static_cast<std::size_t>(a);
          vulnerable_biomass =
              vulnerable_biomass + n[i] * T(weight[i]) * sel[i];
          pred_age_comp[i] = n[i] * sel[i];
          selected_numbers_sum = selected_numbers_sum + pred_age_comp[i];
        }

        const T index_hat = fleet_q * vulnerable_biomass;

        if (obs.index > 0.0) {
          const T z =
              (log_t(T(obs.index)) - log_t(max_t(index_hat, min_positive))) /
              sigma_log_index;
          nll = nll + T(0.5) * square_t(z);
        }

        const T catch_hat = total_catch_hat * T(fleet_catch_share(obs.fleet));
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

      std::array<T, kAges> next{};
      next[0] = r0 * exp_t(rec_dev);

      for (int a = 1; a < kAges; ++a) {
        const auto prev = static_cast<std::size_t>(a - 1);
        const T avg_sel_prev =
            T(0.5) * (sel_longline[prev] + sel_purse_seine[prev]);
        const T f_prev = fbar * avg_sel_prev;
        const T z_prev = m + f_prev;
        next[static_cast<std::size_t>(a)] = n[prev] * exp_t(-z_prev);
      }

      const auto last = static_cast<std::size_t>(kAges - 1);
      const T avg_sel_last =
          T(0.5) * (sel_longline[last] + sel_purse_seine[last]);
      const T f_last = fbar * avg_sel_last;
      const T z_last = m + f_last;
      next[last] = next[last] + n[last] * exp_t(-z_last);

      n = next;
    }

    return nll;
  }

  std::size_t n_years() const { return unique_years().size(); }

  std::vector<int> unique_years() const {
    std::vector<int> years;
    for (const auto &obs : observations_) {
      if (years.empty() || years.back() != obs.year) {
        years.push_back(obs.year);
      }
    }
    return years;
  }

  double fleet_catch_share(const std::string &fleet) const {
    double fleet_total = 0.0;
    double all_total = 0.0;

    for (const auto &obs : observations_) {
      all_total += obs.catch_mt;
      if (obs.fleet == fleet) {
        fleet_total += obs.catch_mt;
      }
    }

    return all_total > 0.0 ? fleet_total / all_total : 0.5;
  }

private:
  std::vector<FleetObservation> observations_;
};

} // namespace pifsc_bigeye_tuna

using pifsc_bigeye_tuna::BigeyeQuadraObjective;
using pifsc_bigeye_tuna::FleetObservation;
