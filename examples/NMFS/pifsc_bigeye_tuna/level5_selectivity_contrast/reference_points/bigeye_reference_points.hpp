#pragma once

#include "../quadra/bigeye_age_structured.hpp"

#include "../../../../../core/optimizer.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

namespace pifsc_bigeye_tuna {

struct BigeyeReferencePoints {
  double r0 = std::numeric_limits<double>::quiet_NaN();
  double fbar = std::numeric_limits<double>::quiet_NaN();
  double q = std::numeric_limits<double>::quiet_NaN();
  double sel_a50 = std::numeric_limits<double>::quiet_NaN();
  double sel_slope = std::numeric_limits<double>::quiet_NaN();

  double spr0 = std::numeric_limits<double>::quiet_NaN();
  double spr_current = std::numeric_limits<double>::quiet_NaN();
  double spr_ratio = std::numeric_limits<double>::quiet_NaN();
  double f40 = std::numeric_limits<double>::quiet_NaN();
  double fbar_over_f40 = std::numeric_limits<double>::quiet_NaN();

  double ssb0_proxy = std::numeric_limits<double>::quiet_NaN();
  double ssb40_proxy = std::numeric_limits<double>::quiet_NaN();
  double terminal_biomass = std::numeric_limits<double>::quiet_NaN();
  double terminal_ssb_proxy = std::numeric_limits<double>::quiet_NaN();
  double terminal_ssb_over_ssb0 = std::numeric_limits<double>::quiet_NaN();
  double terminal_ssb_over_ssb40 = std::numeric_limits<double>::quiet_NaN();
};

inline AgeStructuredParams
bigeye_params_from_fit(const quadra::OptResult &fit) {
  AgeStructuredParams params;

  if (fit.par.size() >= 3) {
    params.log_r0 = fit.par[0];
    params.log_fbar = fit.par[1];
    params.log_q = fit.par[2];
  }

  if (fit.par.size() >= 5) {
    params.sel_a50 = 1.0 + 9.0 / (1.0 + std::exp(-fit.par[3]));
    params.sel_slope = std::exp(fit.par[4]);
  }

  return params;
}

inline double bigeye_spr_proxy(double fbar, double sel_a50, double sel_slope) {
  const double m = 0.18;
  const auto weight = default_weight_at_age();
  const auto maturity = default_maturity_at_age();

  std::vector<double> n(static_cast<std::size_t>(kAges), 0.0);
  n[0] = 1.0;

  for (int a = 1; a < kAges; ++a) {
    const double age_prev = static_cast<double>(a);
    const double sel_prev = logistic_selectivity(age_prev, sel_a50, sel_slope);
    const double z_prev = m + fbar * sel_prev;

    n[static_cast<std::size_t>(a)] =
        n[static_cast<std::size_t>(a - 1)] * std::exp(-z_prev);
  }

  {
    const int last_age = kAges;
    const double sel_last =
        logistic_selectivity(static_cast<double>(last_age), sel_a50, sel_slope);
    const double z_last = m + fbar * sel_last;
    const double denom = std::max(1.0e-12, 1.0 - std::exp(-z_last));

    n[static_cast<std::size_t>(kAges - 1)] =
        n[static_cast<std::size_t>(kAges - 1)] / denom;
  }

  double spr = 0.0;
  for (int a = 0; a < kAges; ++a) {
    const auto i = static_cast<std::size_t>(a);
    spr += n[i] * weight[i] * maturity[i];
  }

  return spr;
}

inline double bigeye_find_f_for_spr_ratio(double target_ratio, double sel_a50,
                                          double sel_slope) {
  const double spr0 = bigeye_spr_proxy(0.0, sel_a50, sel_slope);
  if (!std::isfinite(spr0) || spr0 <= 0.0) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  auto ratio_at = [&](double f) {
    return bigeye_spr_proxy(f, sel_a50, sel_slope) / spr0;
  };

  double lo = 0.0;
  double hi = 0.05;

  while (ratio_at(hi) > target_ratio && hi < 10.0) {
    hi *= 2.0;
  }

  if (hi >= 10.0 && ratio_at(hi) > target_ratio) {
    return std::numeric_limits<double>::quiet_NaN();
  }

  for (int iter = 0; iter < 80; ++iter) {
    const double mid = 0.5 * (lo + hi);
    if (ratio_at(mid) > target_ratio)
      lo = mid;
    else
      hi = mid;
  }

  return 0.5 * (lo + hi);
}

inline BigeyeReferencePoints
compute_bigeye_reference_points(const std::vector<Observation> &observations,
                                const quadra::OptResult &fit) {
  BigeyeReferencePoints out;

  if (fit.par.size() < 5) {
    return out;
  }

  const AgeStructuredParams params = bigeye_params_from_fit(fit);

  out.r0 = std::exp(params.log_r0);
  out.fbar = std::exp(params.log_fbar);
  out.q = std::exp(params.log_q);
  out.sel_a50 = params.sel_a50;
  out.sel_slope = params.sel_slope;

  out.spr0 = bigeye_spr_proxy(0.0, out.sel_a50, out.sel_slope);
  out.spr_current = bigeye_spr_proxy(out.fbar, out.sel_a50, out.sel_slope);

  if (std::isfinite(out.spr0) && out.spr0 > 0.0) {
    out.spr_ratio = out.spr_current / out.spr0;
    out.ssb0_proxy = out.r0 * out.spr0;
    out.ssb40_proxy = 0.4 * out.ssb0_proxy;
  }

  out.f40 = bigeye_find_f_for_spr_ratio(0.4, out.sel_a50, out.sel_slope);
  if (std::isfinite(out.f40) && out.f40 > 0.0) {
    out.fbar_over_f40 = out.fbar / out.f40;
  }

  const auto rows =
      run_deterministic_age_structured_model(observations, params);
  if (!rows.empty()) {
    const auto &terminal = rows.back();
    out.terminal_biomass = terminal.total_biomass;
    out.terminal_ssb_proxy = terminal.ssb_proxy;

    if (std::isfinite(out.ssb0_proxy) && out.ssb0_proxy > 0.0) {
      out.terminal_ssb_over_ssb0 = out.terminal_ssb_proxy / out.ssb0_proxy;
    }

    if (std::isfinite(out.ssb40_proxy) && out.ssb40_proxy > 0.0) {
      out.terminal_ssb_over_ssb40 = out.terminal_ssb_proxy / out.ssb40_proxy;
    }
  }

  return out;
}

inline void write_bigeye_reference_points_csv(const std::string &path,
                                              const BigeyeReferencePoints &rp) {
  std::ofstream out(path);
  out << "quantity,value,note\n";
  out << "R0," << rp.r0 << ",estimated unfished recruitment scale\n";
  out << "Fbar," << rp.fbar
      << ",estimated fully selected fishing mortality proxy\n";
  out << "q," << rp.q << ",catchability estimate\n";
  out << "selectivity_a50," << rp.sel_a50
      << ",estimated logistic selectivity age at 50 percent selected\n";
  out << "selectivity_slope," << rp.sel_slope
      << ",estimated logistic selectivity slope\n";
  out << "SPR0_proxy," << rp.spr0
      << ",unfished spawning output per recruit proxy\n";
  out << "SPR_current_proxy," << rp.spr_current
      << ",spawning output per recruit proxy at estimated Fbar\n";
  out << "SPR_ratio," << rp.spr_ratio
      << ",SPR_current_proxy divided by SPR0_proxy\n";
  out << "F40_proxy," << rp.f40
      << ",Fbar value giving 40 percent SPR proxy under fitted selectivity\n";
  out << "Fbar_over_F40_proxy," << rp.fbar_over_f40
      << ",estimated Fbar relative to F40 proxy\n";
  out << "SSB0_proxy," << rp.ssb0_proxy << ",R0 multiplied by SPR0_proxy\n";
  out << "SSB40_proxy," << rp.ssb40_proxy << ",40 percent of SSB0_proxy\n";
  out << "terminal_biomass," << rp.terminal_biomass
      << ",terminal deterministic total biomass from fitted trajectory\n";
  out << "terminal_SSB_proxy," << rp.terminal_ssb_proxy
      << ",terminal deterministic spawning biomass proxy from fitted "
         "trajectory\n";
  out << "terminal_SSB_over_SSB0_proxy," << rp.terminal_ssb_over_ssb0
      << ",terminal SSB proxy relative to unfished SSB proxy\n";
  out << "terminal_SSB_over_SSB40_proxy," << rp.terminal_ssb_over_ssb40
      << ",terminal SSB proxy relative to SSB40 proxy\n";
}

inline void
write_bigeye_reference_points_csv(const std::string &path,
                                  const std::vector<Observation> &observations,
                                  const quadra::OptResult &fit) {
  write_bigeye_reference_points_csv(
      path, compute_bigeye_reference_points(observations, fit));
}

} // namespace pifsc_bigeye_tuna

using pifsc_bigeye_tuna::BigeyeReferencePoints;
using pifsc_bigeye_tuna::compute_bigeye_reference_points;
using pifsc_bigeye_tuna::write_bigeye_reference_points_csv;
