#pragma once

#include "../quadra/opakapaka_model.hpp"

#include <cmath>
#include <fstream>
#include <limits>
#include <string>
#include <vector>

namespace opakapaka_example {

struct OpakapakaReferencePoints {
  double q = std::numeric_limits<double>::quiet_NaN();
  double r = std::numeric_limits<double>::quiet_NaN();
  double K = std::numeric_limits<double>::quiet_NaN();

  double B_MSY = std::numeric_limits<double>::quiet_NaN();
  double F_MSY = std::numeric_limits<double>::quiet_NaN();
  double MSY = std::numeric_limits<double>::quiet_NaN();

  double B_terminal = std::numeric_limits<double>::quiet_NaN();
  double B_terminal_over_B_MSY = std::numeric_limits<double>::quiet_NaN();
  double F_status_quo = std::numeric_limits<double>::quiet_NaN();
  double F_status_quo_over_F_MSY = std::numeric_limits<double>::quiet_NaN();
};

inline OpakapakaReferencePoints
compute_opakapaka_reference_points(const quadra::OptResult &fit,
                                   const std::vector<Observation> &data) {
  OpakapakaReferencePoints out;

  if (fit.par.size() < 3 || fit.u_hat.empty()) {
    return out;
  }

  out.q = std::exp(fit.par.at(0));
  out.r = std::exp(fit.par.at(1));
  out.K = std::exp(fit.par.at(2));

  out.B_MSY = 0.5 * out.K;
  out.F_MSY = 0.5 * out.r;
  out.MSY = 0.25 * out.r * out.K;

  out.B_terminal = std::exp(fit.u_hat.back());
  if (std::isfinite(out.B_MSY) && out.B_MSY > 0.0) {
    out.B_terminal_over_B_MSY = out.B_terminal / out.B_MSY;
  }

  if (!data.empty() && std::isfinite(out.B_terminal) && out.B_terminal > 0.0) {
    const double recent_catch = data.back().catch_mt;
    out.F_status_quo = recent_catch / out.B_terminal;
    if (std::isfinite(out.F_MSY) && out.F_MSY > 0.0) {
      out.F_status_quo_over_F_MSY = out.F_status_quo / out.F_MSY;
    }
  }

  return out;
}

inline void
write_opakapaka_reference_points_csv(const std::string &path,
                                     const OpakapakaReferencePoints &rp) {
  std::ofstream out(path);
  out << "quantity,value,note\n";
  out << "q," << rp.q << ",catchability estimate\n";
  out << "r," << rp.r << ",intrinsic growth rate estimate\n";
  out << "K," << rp.K << ",carrying capacity estimate\n";
  out << "B_MSY," << rp.B_MSY
      << ",Schaefer surplus-production biomass at MSY equals K/2\n";
  out << "F_MSY," << rp.F_MSY
      << ",Schaefer surplus-production fishing mortality proxy equals r/2\n";
  out << "MSY," << rp.MSY
      << ",Schaefer surplus-production maximum sustainable yield equals "
         "r*K/4\n";
  out << "B_terminal," << rp.B_terminal << ",terminal fitted biomass state\n";
  out << "B_terminal_over_B_MSY," << rp.B_terminal_over_B_MSY
      << ",terminal biomass relative to B_MSY\n";
  out << "F_status_quo," << rp.F_status_quo
      << ",recent catch divided by terminal biomass\n";
  out << "F_status_quo_over_F_MSY," << rp.F_status_quo_over_F_MSY
      << ",status quo fishing mortality proxy relative to F_MSY\n";
}

inline void
write_opakapaka_reference_points_csv(const std::string &path,
                                     const quadra::OptResult &fit,
                                     const std::vector<Observation> &data) {
  write_opakapaka_reference_points_csv(
      path, compute_opakapaka_reference_points(fit, data));
}

} // namespace opakapaka_example

using opakapaka_example::compute_opakapaka_reference_points;
using opakapaka_example::OpakapakaReferencePoints;
using opakapaka_example::write_opakapaka_reference_points_csv;
