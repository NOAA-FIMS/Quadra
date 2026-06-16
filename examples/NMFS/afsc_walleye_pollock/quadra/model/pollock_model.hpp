#pragma once

#include "pollock_constants.hpp"

#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>

struct Obs
{
  int year;
  double catch_mt;
  double index;
  std::vector<double> age;
};

struct PollockModel
{
  explicit PollockModel(std::vector<Obs> obs) : obs_(std::move(obs)) {}

  template <class AD>
  AD operator()(const std::vector<AD> &p) const
  {
    const AD log_r0 = p[0];
    const AD log_fbar = p[1];

    const AD r0 = exp(log_r0);
    const AD fbar = exp(log_fbar);

    // Assessment-like initialization:
    // derive initial numbers-at-age from the same unfished recruitment scale
    // that drives the stock-recruit curve. This removes the artificial
    // R0/N0 conflict from the synthetic scaffold.
    const AD n0 = r0;

    // Hold catchability fixed in this scaffold-level scaling experiment.
    // This isolates recruitment/Laplace behavior from q-abundance
    // confounding after selectivity has already been fixed.
    const AD q = exp(AD(-8.78));

    // Hold selectivity fixed in this scaffold-level scaling experiment.
    // This isolates recruitment/Laplace behavior from selectivity-q-abundance
    // confounding.
    const AD sel_a50 = AD(4.0);
    const AD sel_slope = AD(1.0);

    constexpr int A = 7;
    const double weight[A] = {0.20, 0.45, 0.75, 1.10, 1.45, 1.75, 2.00};
    const double maturity[A] = {0.00, 0.10, 0.45, 0.80, 0.95, 1.00, 1.00};
    const double M = 0.25;

    AD nll = AD(0.0);
    nll += AD(0.5) * pow((log_r0 - AD(8.0)) / AD(4.0), 2.0);
    nll += AD(0.5) * pow((log_fbar - AD(-3.7)) / AD(3.0), 2.0);

    // Equilibrium numbers-at-age from the R0 recruitment scale.
    // Ages 1..A-1 follow survivorship; the terminal age is a plus group
    // accumulating survivors from all older cohorts.
    std::vector<AD> N(pollock::n_ages);
    const AD surv = exp(-AD(pollock::natural_mortality));
    N[0] = r0;
    for (int a = 1; a < pollock::n_ages - 1; ++a)
    {
      N[a] = N[a - 1] * surv;
    }
    N[pollock::n_ages - 1] = N[pollock::n_ages - 2] * surv / (AD(1.0) - surv);

    const AD rec_sigma = AD(0.15);
    const AD rec_rho = AD(0.60);
    const AD rec_stationary_sigma =
        rec_sigma / sqrt(AD(1.0) - rec_rho * rec_rho);
    const AD index_sigma = AD(0.30);
    const AD catch_sigma = AD(0.25);
#ifdef WALLEYE_POLLOCK_FIT_CATCH_LIKELIHOOD
    const AD catch_w = AD(1.0);
#else
    const AD catch_w = AD(0.0);
#endif
    const AD age_w = AD(0.0);

    for (std::size_t y = 0; y < obs_.size(); ++y)
    {
      const std::size_t rec_offset = 2;
      const bool has_rec_dev = (p.size() > rec_offset + y);
      const AD rec_dev =
          has_rec_dev ? p[rec_offset + y] : AD(0.0);

      if (has_rec_dev)
      {
        if (y == 0)
        {
          // Stationary AR(1) prior for the initial recruitment deviation.
          nll += AD(0.5) * pow(rec_dev / rec_stationary_sigma, 2.0) +
                 log(rec_stationary_sigma);
        }
        else
        {
          const bool has_prev_rec_dev = (p.size() > rec_offset + y - 1);
          const AD prev_rec_dev =
              has_prev_rec_dev ? p[rec_offset + y - 1] : AD(0.0);
          const AD innovation = rec_dev - rec_rho * prev_rec_dev;

          nll += AD(0.5) * pow(innovation / rec_sigma, 2.0) +
                 log(rec_sigma);
        }
      }

      AD biomass = AD(0.0);
      AD ssb = AD(0.0);
      AD pred_catch = AD(0.0);
      std::vector<AD> caa(pollock::n_ages);

      for (int a = 0; a < pollock::n_ages; ++a)
      {
        const AD sel = AD(1.0) / (AD(1.0) + exp(-sel_slope * (AD(a + 1) - sel_a50)));
        const AD Z = AD(pollock::natural_mortality) + fbar * sel;
        biomass += N[a] * AD(pollock::weight_at_age[a]);
        ssb += N[a] * AD(pollock::weight_at_age[a] * pollock::maturity_at_age[a]);
        caa[a] = N[a] * (fbar * sel / Z) * (AD(1.0) - exp(-Z)) * AD(pollock::weight_at_age[a]);
        pred_catch += caa[a];
      }

      const AD pred_index = q * biomass;

      nll += AD(0.5) * pow((log(AD(obs_[y].index) + AD(1e-12)) -
                            log(pred_index + AD(1e-12))) /
                               index_sigma,
                           2.0);
      if (catch_w > AD(0.0))
      {
        nll += catch_w *
               AD(0.5) *
               pow((log(AD(obs_[y].catch_mt) + AD(1e-12)) -
                    log(pred_catch + AD(1e-12))) /
                       catch_sigma,
                   2.0);
      }

      if (age_w > AD(0.0))
      {
        for (int a = 0; a < pollock::n_ages; ++a)
        {
          const AD pred_p = caa[a] / (pred_catch + AD(1e-12));
          nll -= age_w * AD(obs_[y].age[a]) * log(pred_p + AD(1e-12));
        }
      }

      std::vector<AD> next(pollock::n_ages);
    // Treat observed catch as the removals driver for this synthetic
    // assessment scaffold. fbar still controls age-specific selectivity and
    // relative exploitation, but total removals are scaled toward observed
    // catch rather than forcing a single constant F to fit the catch series.
    const AD catch_scale_raw =
        AD(obs_[y].catch_mt) / (pred_catch + AD(1.0e-12));
    const AD catch_scale =
        (catch_scale_raw < AD(0.95)) ? catch_scale_raw : AD(0.95);

    for (int a = 0; a < pollock::n_ages; ++a)
    {
      const AD catch_number =
          catch_scale * caa[a] / (AD(pollock::weight_at_age[a]) + AD(1.0e-12));
      N[a] = (N[a] > catch_number) ? (N[a] - catch_number) : AD(1.0e-12);
    }

    // Ricker-style stock-recruitment relationship.
    //
    // This synthetic scaffold anchors the curve so that R(B0) = R0 at an
    // approximate unfished spawning biomass B0. Recruitment deviations remain
    // multiplicative lognormal random effects around the stock-recruit curve.
    const AD b0 = r0 * AD(4.0);
    const AD beta = AD(1.0) / (b0 + AD(1.0e-12));
    const AD alpha = r0 * exp(beta * b0) / (b0 + AD(1.0e-12));
    const AD recruitment =
        alpha * ssb * exp(-beta * ssb + rec_dev);

    next[0] = recruitment;
      for (int a = 1; a < pollock::n_ages; ++a)
        next[a] = N[a - 1] * exp(-AD(pollock::natural_mortality));
      next[pollock::n_ages - 1] += N[pollock::n_ages - 1] * exp(-AD(pollock::natural_mortality));
      N = next;
    }

    return nll;
  }

  std::vector<Obs> obs_;
};
