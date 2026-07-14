#pragma once

#include "../../../../../core/optimizer.hpp"
#include "../quadra/bigeye_age_structured.hpp"

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <numeric>
#include <stdexcept>
#include <string>

namespace pifsc_bigeye_tuna {

inline void
write_level13_initial_numbers_diagnostics(const std::string &text_path,
                                          const std::string &csv_path,
                                          const quadra::OptResult &fit) {
  constexpr int kBaseFixed = 3;
  constexpr int kInitialOffset = kBaseFixed;

  if (fit.par.size() < static_cast<std::size_t>(kBaseFixed + kAges))
    throw std::runtime_error(
        "Level 20 initial-number diagnostics expected base fixed effects plus "
        "age-specific initial deviations");

  const double log_r0 = fit.par[0];
  const double log_m = std::log(0.45);
  const double r0 = std::exp(log_r0);
  const double m = std::exp(log_m);

  const auto weight = default_weight_at_age();
  const auto maturity = default_maturity_at_age();

  std::array<double, kAges> equilibrium_n{};
  equilibrium_n[0] = r0;
  for (int a = 1; a < kAges; ++a)
    equilibrium_n[static_cast<std::size_t>(a)] =
        equilibrium_n[static_cast<std::size_t>(a - 1)] * std::exp(-m);

  // Plus group.
  equilibrium_n[static_cast<std::size_t>(kAges - 1)] =
      equilibrium_n[static_cast<std::size_t>(kAges - 1)] / (1.0 - std::exp(-m));

  std::array<double, kAges> fitted_n{};
  std::array<double, kAges> multipliers{};
  for (int a = 0; a < kAges; ++a) {
    const auto i = static_cast<std::size_t>(a);
    multipliers[i] =
        std::exp(fit.par[static_cast<std::size_t>(kInitialOffset + a)]);
    fitted_n[i] = equilibrium_n[i] * multipliers[i];
  }

  auto biomass = [&](const std::array<double, kAges> &n) {
    double out = 0.0;
    for (int a = 0; a < kAges; ++a)
      out +=
          n[static_cast<std::size_t>(a)] * weight[static_cast<std::size_t>(a)];
    return out;
  };

  auto ssb = [&](const std::array<double, kAges> &n) {
    double out = 0.0;
    for (int a = 0; a < kAges; ++a)
      out += n[static_cast<std::size_t>(a)] *
             weight[static_cast<std::size_t>(a)] *
             maturity[static_cast<std::size_t>(a)];
    return out;
  };

  const double eq_total_n =
      std::accumulate(equilibrium_n.begin(), equilibrium_n.end(), 0.0);
  const double fit_total_n =
      std::accumulate(fitted_n.begin(), fitted_n.end(), 0.0);

  const double eq_biomass = biomass(equilibrium_n);
  const double fit_biomass = biomass(fitted_n);
  const double eq_ssb = ssb(equilibrium_n);
  const double fit_ssb = ssb(fitted_n);

  const int plus = kAges - 1;
  const double eq_plus_n = equilibrium_n[static_cast<std::size_t>(plus)];
  const double fit_plus_n = fitted_n[static_cast<std::size_t>(plus)];
  const double eq_plus_biomass =
      eq_plus_n * weight[static_cast<std::size_t>(plus)];
  const double fit_plus_biomass =
      fit_plus_n * weight[static_cast<std::size_t>(plus)];

  std::ofstream txt(text_path);
  if (!txt)
    throw std::runtime_error(
        "Could not open initial numbers diagnostic text: " + text_path);

  std::ofstream csv(csv_path);
  if (!csv)
    throw std::runtime_error("Could not open initial numbers diagnostic CSV: " +
                             csv_path);

  txt << std::setprecision(15);
  csv << std::setprecision(15);

  txt << "Level 20 Initial Numbers Diagnostics\n";
  txt << "====================================\n\n";
  txt << "Interpretation\n";
  txt << "--------------\n";
  txt << "Equilibrium initial numbers are the numbers-at-age implied by R0 and "
         "fixed M\n";
  txt << "before estimating the age-specific initial-number deviations. Fitted "
         "initial\n";
  txt << "numbers multiply that equilibrium vector by "
         "exp(init_log_number_dev_age_a).\n";
  txt << "Large reductions in old ages or the plus group indicate that "
         "recruitment\n";
  txt << "deviations were likely compensating for an incorrect initial age "
         "structure.\n\n";

  txt << "Summary\n";
  txt << "-------\n";
  txt << "r0:                         " << r0 << "\n";
  txt << "m_fixed:                    " << m << "\n";
  txt << "equilibrium_total_n:         " << eq_total_n << "\n";
  txt << "fitted_total_n:              " << fit_total_n << "\n";
  txt << "fitted_over_equilibrium_n:   " << fit_total_n / eq_total_n << "\n";
  txt << "equilibrium_biomass:         " << eq_biomass << "\n";
  txt << "fitted_biomass:              " << fit_biomass << "\n";
  txt << "fitted_over_equilibrium_bio: " << fit_biomass / eq_biomass << "\n";
  txt << "equilibrium_ssb_proxy:       " << eq_ssb << "\n";
  txt << "fitted_ssb_proxy:            " << fit_ssb << "\n";
  txt << "fitted_over_equilibrium_ssb: " << fit_ssb / eq_ssb << "\n";
  txt << "equilibrium_plus_n_share:    " << eq_plus_n / eq_total_n << "\n";
  txt << "fitted_plus_n_share:         " << fit_plus_n / fit_total_n << "\n";
  txt << "equilibrium_plus_bio_share:  " << eq_plus_biomass / eq_biomass
      << "\n";
  txt << "fitted_plus_bio_share:       " << fit_plus_biomass / fit_biomass
      << "\n\n";

  txt << "Rows\n";
  txt << "----\n";
  txt << "age,equilibrium_n,init_log_dev,multiplier,fitted_n,"
         "equilibrium_biomass,fitted_biomass,equilibrium_n_share,"
         "fitted_n_share,equilibrium_biomass_share,fitted_biomass_share,"
         "weight,maturity\n";

  csv << "section,metric,target,value,extra\n";
  csv << "summary,r0,," << r0 << ",\n";
  csv << "summary,m_fixed,," << m << ",\n";
  csv << "summary,equilibrium_total_n,," << eq_total_n << ",\n";
  csv << "summary,fitted_total_n,," << fit_total_n << ",\n";
  csv << "summary,fitted_over_equilibrium_n,," << fit_total_n / eq_total_n
      << ",\n";
  csv << "summary,equilibrium_biomass,," << eq_biomass << ",\n";
  csv << "summary,fitted_biomass,," << fit_biomass << ",\n";
  csv << "summary,fitted_over_equilibrium_bio,," << fit_biomass / eq_biomass
      << ",\n";
  csv << "summary,equilibrium_ssb_proxy,," << eq_ssb << ",\n";
  csv << "summary,fitted_ssb_proxy,," << fit_ssb << ",\n";
  csv << "summary,fitted_over_equilibrium_ssb,," << fit_ssb / eq_ssb << ",\n";
  csv << "summary,equilibrium_plus_n_share,," << eq_plus_n / eq_total_n
      << ",\n";
  csv << "summary,fitted_plus_n_share,," << fit_plus_n / fit_total_n << ",\n";
  csv << "summary,equilibrium_plus_bio_share,," << eq_plus_biomass / eq_biomass
      << ",\n";
  csv << "summary,fitted_plus_bio_share,," << fit_plus_biomass / fit_biomass
      << ",\n";

  for (int a = 0; a < kAges; ++a) {
    const auto i = static_cast<std::size_t>(a);
    const double init_log_dev =
        fit.par[static_cast<std::size_t>(kInitialOffset + a)];
    const double eq_bio_a = equilibrium_n[i] * weight[i];
    const double fit_bio_a = fitted_n[i] * weight[i];

    txt << (a + 1) << "," << equilibrium_n[i] << "," << init_log_dev << ","
        << multipliers[i] << "," << fitted_n[i] << "," << eq_bio_a << ","
        << fit_bio_a << "," << equilibrium_n[i] / eq_total_n << ","
        << fitted_n[i] / fit_total_n << "," << eq_bio_a / eq_biomass << ","
        << fit_bio_a / fit_biomass << "," << weight[i] << "," << maturity[i]
        << "\n";

    csv << "age,equilibrium_n,age_" << (a + 1) << "," << equilibrium_n[i]
        << ",\n";
    csv << "age,init_log_dev,age_" << (a + 1) << "," << init_log_dev << ",\n";
    csv << "age,multiplier,age_" << (a + 1) << "," << multipliers[i] << ",\n";
    csv << "age,fitted_n,age_" << (a + 1) << "," << fitted_n[i] << ",\n";
    csv << "age,equilibrium_biomass,age_" << (a + 1) << "," << eq_bio_a
        << ",\n";
    csv << "age,fitted_biomass,age_" << (a + 1) << "," << fit_bio_a << ",\n";
    csv << "age,equilibrium_n_share,age_" << (a + 1) << ","
        << equilibrium_n[i] / eq_total_n << ",\n";
    csv << "age,fitted_n_share,age_" << (a + 1) << ","
        << fitted_n[i] / fit_total_n << ",\n";
    csv << "age,equilibrium_biomass_share,age_" << (a + 1) << ","
        << eq_bio_a / eq_biomass << ",\n";
    csv << "age,fitted_biomass_share,age_" << (a + 1) << ","
        << fit_bio_a / fit_biomass << ",\n";
    csv << "age,weight,age_" << (a + 1) << "," << weight[i] << ",\n";
    csv << "age,maturity,age_" << (a + 1) << "," << maturity[i] << ",\n";
  }
}

} // namespace pifsc_bigeye_tuna

using pifsc_bigeye_tuna::write_level13_initial_numbers_diagnostics;
