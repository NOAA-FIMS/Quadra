#pragma once
#include "bigeye_level21_m_at_age_helpers.hpp"

#include "bigeye_level21_parameter_layout.hpp"

#include "../objective/bigeye_quadra_objective.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace pifsc_bigeye_tuna {

using namespace level21_layout;

inline void write_level16_purse_seine_prediction_decomposition(
    const std::string &txt_path, const std::string &csv_path,
    const BigeyeQuadraObjective &objective, const quadra::OptResult &fit) {
  std::ofstream txt(txt_path);
  std::ofstream csv(csv_path);
  if (!txt || !csv) {
    throw std::runtime_error(
        "Cannot write purse-seine prediction decomposition");
  }

  txt << std::setprecision(15);
  csv << std::setprecision(15);
  const double min_positive = 1.0e-12;

  const double r0 = std::exp(fit.par[0]);

  const auto weight = default_weight_at_age();
  const auto maturity = default_maturity_at_age();
  const double fbar = std::exp(fit.par[1]);
  const double q_purse_seine = std::exp(fit.par[2]);
  const double adult_m = 0.45;
  const auto m_at_age = level21_m_helpers::m_at_age_from_level21_par(fit.par);
  std::array<double, kAges> sel_longline{};
  std::array<double, kAges> sel_purse_seine{};
  for (int a = 0; a < kAges; ++a) {
    sel_longline[static_cast<std::size_t>(a)] =
        1.0 /
        (1.0 + std::exp(-fit.par[level21_m_helpers::kLonglineSelOffset + a]));
    sel_purse_seine[static_cast<std::size_t>(a)] =
        1.0 /
        (1.0 + std::exp(-fit.par[level21_m_helpers::kPurseSeineSelOffset + a]));
  }

  std::array<double, kAges> n{};
  n[0] = r0;
  for (int a = 1; a < kAges; ++a) {
    n[static_cast<std::size_t>(a)] =
        n[static_cast<std::size_t>(a - 1)] *
        std::exp(-m_at_age[static_cast<std::size_t>(a - 1)]);
  }
  n[static_cast<std::size_t>(kAges - 1)] =
      n[static_cast<std::size_t>(kAges - 1)] /
      (1.0 - std::exp(-m_at_age[static_cast<std::size_t>(kAges - 1)]));

  for (int a = 0; a < kAges; ++a) {
    n[static_cast<std::size_t>(a)] *=
        std::exp(fit.par[level21_m_helpers::kInitialDevOffset + a]);
  }

  txt << "Level 21 Purse-Seine Prediction Decomposition\n";
  txt << "=============================================\n\n";
  txt << "Purpose\n";
  txt << "-------\n";
  txt << "Decompose purse-seine predicted age composition into "
         "numbers-at-age,\n";
  txt << "purse-seine selectivity, selected numbers, predicted composition,\n";
  txt << "observed composition, and residual.\n\n";

  txt << "Fixed quantities\n";
  txt << "----------------\n";
  txt << "r0: " << r0 << "\n";
  txt << "m_young: " << m_at_age[0] << "\n";
  txt << "m_adult: " << m_at_age[3] << "\n";
  txt << "m_old: " << m_at_age[kAges - 1] << "\n";
  txt << "fbar: " << fbar << "\n";
  txt << "q_purse_seine: " << q_purse_seine << "\n\n";

  csv << "section,year,age,n_at_age,purse_seine_selectivity,selected_numbers,"
         "predicted_comp,observed_comp,residual,abs_residual,selected_share,"
         "fleet,index,catch_mt\n";

  txt << "Rows\n";
  txt << "----\n";
  txt << "year,age,n_at_age,purse_seine_selectivity,selected_numbers,"
         "predicted_comp,observed_comp,residual,abs_residual,selected_share\n";

  struct ResidualRow {
    int year = 0;
    int age = 0;
    double observed = 0.0;
    double predicted = 0.0;
    double residual = 0.0;
    double abs_residual = 0.0;
    double n_at_age = 0.0;
    double sel = 0.0;
    double selected = 0.0;
  };

  std::vector<ResidualRow> residual_rows;

  for (std::size_t t = 0; t < objective.fleet_observations().size(); ++t) {
    const auto &obs = objective.fleet_observations()[t];

    if (obs.fleet == "purse_seine") {
      std::array<double, kAges> selected{};
      double selected_sum = 0.0;

      for (int a = 0; a < kAges; ++a) {
        const std::size_t i = static_cast<std::size_t>(a);
        selected[i] = n[i] * sel_purse_seine[i];
        selected_sum += selected[i];
      }

      for (int a = 0; a < kAges; ++a) {
        const std::size_t i = static_cast<std::size_t>(a);
        const double pred = selected[i] / std::max(selected_sum, min_positive);
        const double obs_a = obs.age_comp[i];
        const double resid = obs_a - pred;
        const double abs_resid = std::abs(resid);

        txt << obs.year << "," << (a + 1) << "," << n[i] << ","
            << sel_purse_seine[i] << "," << selected[i] << "," << pred << ","
            << obs_a << "," << resid << "," << abs_resid << "," << pred << "\n";

        csv << "purse_seine_decomposition," << obs.year << "," << (a + 1) << ","
            << n[i] << "," << sel_purse_seine[i] << "," << selected[i] << ","
            << pred << "," << obs_a << "," << resid << "," << abs_resid << ","
            << pred << "," << obs.fleet << "," << obs.index << ","
            << obs.catch_mt << "\n";

        residual_rows.push_back({obs.year, a + 1, obs_a, pred, resid, abs_resid,
                                 n[i], sel_purse_seine[i], selected[i]});
      }
    }

    std::array<double, kAges> total_sel{};
    for (int a = 0; a < kAges; ++a) {
      const std::size_t i = static_cast<std::size_t>(a);
      total_sel[i] = sel_longline[i] + sel_purse_seine[i];
    }

    std::array<double, kAges> next{};
    const double phi0 =
        level21_bh_sync::spawning_biomass_from_numbers(n, weight, maturity) /
        std::max(r0, level21_bh_sync::bh_min_positive());

    const double spawning_biomass =
        level21_bh_sync::spawning_biomass_from_numbers(n, weight, maturity);
    const double expected_recruitment =
        level21_bh_sync::beverton_holt_recruitment(spawning_biomass, r0, 0.75,
                                                   phi0);
    next[0] = expected_recruitment *
              std::exp(fit.u_hat[static_cast<Eigen::Index>(t)]);
    for (int a = 1; a < kAges; ++a) {
      const std::size_t prev = static_cast<std::size_t>(a - 1);
      const double z_prev = m_at_age[prev] + fbar * total_sel[prev];
      next[static_cast<std::size_t>(a)] = n[prev] * std::exp(-z_prev);
    }
    const std::size_t last = static_cast<std::size_t>(kAges - 1);
    const double z_last = m_at_age[last] + fbar * total_sel[last];
    next[last] += n[last] * std::exp(-z_last);
    n = next;
  }

  std::sort(residual_rows.begin(), residual_rows.end(),
            [](const ResidualRow &a, const ResidualRow &b) {
              return a.abs_residual > b.abs_residual;
            });

  txt << "\nTop residuals\n";
  txt << "-------------\n";
  txt << "year,age,observed,predicted,residual,abs_residual,n_at_age,"
         "purse_seine_selectivity,selected_numbers\n";

  for (std::size_t i = 0; i < std::min<std::size_t>(20, residual_rows.size());
       ++i) {
    const auto &r = residual_rows[i];
    txt << r.year << "," << r.age << "," << r.observed << "," << r.predicted
        << "," << r.residual << "," << r.abs_residual << "," << r.n_at_age
        << "," << r.sel << "," << r.selected << "\n";
  }

  std::array<double, kAges> mean_obs{};
  std::array<double, kAges> mean_pred{};
  std::array<double, kAges> mean_abs{};
  std::array<int, kAges> counts{};

  for (const auto &r : residual_rows) {
    const std::size_t i = static_cast<std::size_t>(r.age - 1);
    mean_obs[i] += r.observed;
    mean_pred[i] += r.predicted;
    mean_abs[i] += r.abs_residual;
    counts[i] += 1;
  }

  txt << "\nMean by age\n";
  txt << "-----------\n";
  txt << "age,mean_observed,mean_predicted,mean_abs_residual,n\n";
  for (int a = 0; a < kAges; ++a) {
    const std::size_t i = static_cast<std::size_t>(a);
    if (counts[i] > 0) {
      mean_obs[i] /= static_cast<double>(counts[i]);
      mean_pred[i] /= static_cast<double>(counts[i]);
      mean_abs[i] /= static_cast<double>(counts[i]);
    }
    txt << (a + 1) << "," << mean_obs[i] << "," << mean_pred[i] << ","
        << mean_abs[i] << "," << counts[i] << "\n";
  }
}

} // namespace pifsc_bigeye_tuna

using pifsc_bigeye_tuna::write_level16_purse_seine_prediction_decomposition;
