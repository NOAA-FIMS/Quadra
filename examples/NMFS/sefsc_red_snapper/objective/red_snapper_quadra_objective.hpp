#pragma once

#include "../quadra/red_snapper_age_structured.hpp"

#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace sefsc_red_snapper {

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
    if (obs > 0.0) {
      nll = nll - T(effective_n * obs) * log_t(max_t(predicted[i], floor));
    }
  }
  return nll;
}

class RedSnapperQuadraObjective {
public:
  explicit RedSnapperQuadraObjective(std::vector<Observation> observations)
      : observations_(std::move(observations)) {}

  template <class T> T operator()(const std::vector<T> &par) const {
    if (par.size() < 5 + observations_.size()) {
      throw std::runtime_error("RedSnapperQuadraObjective expected parameters: "
                               "log_r0, log_fbar, log_q");
    }

    const T log_r0 = par[0];
    const T log_fbar = par[1];
    const T log_q = par[2];
    const T logit_sel_a50 = par[3];
    const T log_sel_slope = par[4];

    const T r0 = exp_t(log_r0);
    const T m = T(0.18);
    const T fbar = exp_t(log_fbar);
    const T q = exp_t(log_q);
    const T sel_a50 = T(1.0) + T(9.0) * invlogit_t(logit_sel_a50);
    const T sel_slope = exp_t(log_sel_slope);

    const T sigma_log_index = T(0.20);
    const T sigma_log_catch = T(0.15);

    const T sigma_rec_dev = T(0.35);
    const double age_comp_effective_n = 2.0;
    const double min_positive = 1.0e-12;

    const auto weight = default_weight_at_age();
    const auto maturity = default_maturity_at_age();

    std::array<T, kAges> selectivity{};
    for (int a = 0; a < kAges; ++a) {
      selectivity[static_cast<std::size_t>(a)] =
          logistic_selectivity_t(T(a + 1), sel_a50, sel_slope);
    }

    std::array<T, kAges> n{};
    n[0] = r0;
    for (int a = 1; a < kAges; ++a) {
      n[static_cast<std::size_t>(a)] =
          n[static_cast<std::size_t>(a - 1)] * exp_t(-m);
    }
    n[static_cast<std::size_t>(kAges - 1)] =
        n[static_cast<std::size_t>(kAges - 1)] / (T(1.0) - exp_t(-m));

    T nll = T(0.0);
    T fixed_prior_nll = T(0.0);
    T rec_prior_nll = T(0.0);
    T index_nll = T(0.0);
    T catch_nll = T(0.0);
    T age_comp_nll_total = T(0.0);

    auto normal_prior = [](const T &x, double mean, double sd) {
      const T z = (x - T(mean)) / T(sd);
      return T(0.5) * z * z;
    };

    fixed_prior_nll =
        fixed_prior_nll + normal_prior(log_r0, std::log(1200.0), 1.0);
    fixed_prior_nll =
        fixed_prior_nll + normal_prior(log_fbar, std::log(0.025), 0.75);
    fixed_prior_nll =
        fixed_prior_nll + normal_prior(log_q, std::log(0.00005), 1.0);
    fixed_prior_nll = fixed_prior_nll + normal_prior(sel_a50, 4.0, 0.75);
    fixed_prior_nll =
        fixed_prior_nll + normal_prior(log_sel_slope, std::log(1.2), 0.35);

    nll = nll + fixed_prior_nll;

    for (std::size_t t = 0; t < observations_.size(); ++t) {

      const auto &obs = observations_[t];

      const T rec_dev = par[5 + t];

      {
        T term = T(0.5) * square_t(rec_dev / sigma_rec_dev);
        rec_prior_nll = rec_prior_nll + term;
        nll = nll + term;
      }
      T biomass = T(0.0);
      for (int a = 0; a < kAges; ++a) {
        biomass = biomass + n[static_cast<std::size_t>(a)] *
                                T(weight[static_cast<std::size_t>(a)]);
      }

      T catch_hat = T(0.0);
      for (int a = 0; a < kAges; ++a) {
        const auto i = static_cast<std::size_t>(a);
        const T f_a = fbar * selectivity[i];
        const T z_a = m + f_a;
        const T harvest_rate = (f_a / z_a) * (T(1.0) - exp_t(-z_a));
        catch_hat = catch_hat + n[i] * T(weight[i]) * harvest_rate;
      }

      const T index_hat = q * biomass;

      if (obs.index > 0.0) {
        const T z =
            (log_t(T(obs.index)) - log_t(max_t(index_hat, min_positive))) /
            sigma_log_index;
        {
          T term = T(0.5) * square_t(z);
          index_nll = index_nll + term;
          nll = nll + term;
        }
      }

      if (obs.catch_mt > 0.0) {
        const T z =
            (log_t(T(obs.catch_mt)) - log_t(max_t(catch_hat, min_positive))) /
            sigma_log_catch;
        {
          T term = T(0.5) * square_t(z);
          catch_nll = catch_nll + term;
          nll = nll + term;
        }
      }

      std::array<T, kAges> pred_age_comp{};
      T selected_numbers_sum = T(0.0);
      for (int a = 0; a < kAges; ++a) {
        const auto i = static_cast<std::size_t>(a);
        pred_age_comp[i] = n[i] * selectivity[i];
        selected_numbers_sum = selected_numbers_sum + pred_age_comp[i];
      }
      for (int a = 0; a < kAges; ++a) {
        const auto i = static_cast<std::size_t>(a);
        pred_age_comp[i] =
            pred_age_comp[i] / max_t(selected_numbers_sum, min_positive);
      }

      {
        T term = age_comp_nll(obs.age_comp, pred_age_comp, age_comp_effective_n,
                              min_positive);
        age_comp_nll_total = age_comp_nll_total + term;
        nll = nll + term;
      }

      std::array<T, kAges> next{};
      next[0] = r0 * exp_t(rec_dev);

      for (int a = 1; a < kAges; ++a) {
        const auto prev = static_cast<std::size_t>(a - 1);
        const T f_prev = fbar * selectivity[prev];
        const T z_prev = m + f_prev;
        next[static_cast<std::size_t>(a)] = n[prev] * exp_t(-z_prev);
      }

      const auto last = static_cast<std::size_t>(kAges - 1);
      const T f_last = fbar * selectivity[last];
      const T z_last = m + f_last;
      next[last] = next[last] + n[last] * exp_t(-z_last);

      n = next;
    }

    return nll;
  }

private:
  std::vector<Observation> observations_;
};


}  // namespace sefsc_red_snapper

using sefsc_red_snapper::RedSnapperQuadraObjective;
