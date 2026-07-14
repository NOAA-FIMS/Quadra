#pragma once

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

inline void write_level18_longline_prediction_decomposition(
    const std::string &txt_path, const std::string &csv_path,
    const BigeyeQuadraObjective &objective, const quadra::OptResult &fit) {
  std::ofstream txt(txt_path);
  std::ofstream csv(csv_path);
  if (!txt || !csv) {
    throw std::runtime_error("Cannot write longline prediction decomposition");
  }

  txt << std::setprecision(15);
  csv << std::setprecision(15);

  constexpr int kBaseFixed = 3;
  constexpr int kLonglineSelOffset = kBaseFixed;
  constexpr int kLonglineSelDevs = kAges;
  constexpr int kInitialDevs = kAges;
  constexpr int kInitialDevOffset = kLonglineSelOffset + kLonglineSelDevs;
  constexpr int kPurseSeineSelOffset = kInitialDevOffset + kInitialDevs;

  const double min_positive = 1.0e-12;
  const double m = 0.45;

  const double r0 = std::exp(fit.par[0]);
  const auto weight = default_weight_at_age();
  const auto maturity = default_maturity_at_age();
  const double steepness = 0.75;
  const double fbar = std::exp(fit.par[1]);
  const double q_purse_seine = std::exp(fit.par[2]);
  const double q_longline = 0.00005;
  std::array<double, kAges> sel_longline{};
  std::array<double, kAges> sel_purse_seine{};
  for (int a = 0; a < kAges; ++a) {
    sel_longline[static_cast<std::size_t>(a)] =
        1.0 / (1.0 + std::exp(-fit.par[kLonglineSelOffset + a]));
    sel_purse_seine[static_cast<std::size_t>(a)] =
        1.0 / (1.0 + std::exp(-fit.par[kPurseSeineSelOffset + a]));
  }

  std::array<double, kAges> n{};
  n[0] = r0;
  for (int a = 1; a < kAges; ++a) {
    n[static_cast<std::size_t>(a)] =
        n[static_cast<std::size_t>(a - 1)] * std::exp(-m);
  }
  n[static_cast<std::size_t>(kAges - 1)] =
      n[static_cast<std::size_t>(kAges - 1)] / (1.0 - std::exp(-m));

  for (int a = 0; a < kAges; ++a) {
    n[static_cast<std::size_t>(a)] *= std::exp(fit.par[kInitialDevOffset + a]);
  }

  txt << "Level 23 Longline Prediction Decomposition\n";
  txt << "==========================================\n\n";
  txt << "Purpose\n";
  txt << "-------\n";
  txt << "Decompose longline predicted age composition into numbers-at-age,\n";
  txt << "longline selectivity, selected numbers, predicted composition,\n";
  txt << "observed composition, and residual.\n\n";

  txt << "Fixed quantities\n";
  txt << "----------------\n";
  txt << "r0: " << r0 << "\n";
  txt << "m: " << m << "\n";
  txt << "fbar: " << fbar << "\n";
  txt << "q_longline: " << q_longline << "\n";
  txt << "q_purse_seine: " << q_purse_seine << "\n";
  txt << "\n";

  csv << "section,year,age,n_at_age,longline_selectivity,selected_numbers,"
         "predicted_comp,observed_comp,residual,abs_residual,selected_share,"
         "fleet,index,catch_mt\n";

  txt << "Rows\n";
  txt << "----\n";
  txt << "year,age,n_at_age,longline_selectivity,selected_numbers,"
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

  const auto years = objective.unique_years();
  const auto &observations = objective.fleet_observations();
  const double unfished_spawning_biomass =
      level23_bh_sync::spawning_biomass_from_numbers(n, weight, maturity);
  const double phi0 = unfished_spawning_biomass /
                      std::max(r0, level23_bh_sync::bh_min_positive());

  for (std::size_t t = 0; t < years.size(); ++t) {
    for (const auto &obs : observations) {
      if (obs.year != years[t] || obs.fleet != "longline") {
        continue;
      }

      std::array<double, kAges> selected{};
      double selected_sum = 0.0;

      for (int a = 0; a < kAges; ++a) {
        const std::size_t i = static_cast<std::size_t>(a);
        selected[i] = n[i] * sel_longline[i];
        selected_sum += selected[i];
      }

      for (int a = 0; a < kAges; ++a) {
        const std::size_t i = static_cast<std::size_t>(a);
        const double pred = selected[i] / std::max(selected_sum, min_positive);
        const double obs_a = obs.age_comp[i];
        const double resid = obs_a - pred;
        const double abs_resid = std::abs(resid);

        txt << obs.year << "," << (a + 1) << "," << n[i] << ","
            << sel_longline[i] << "," << selected[i] << "," << pred << ","
            << obs_a << "," << resid << "," << abs_resid << "," << pred << "\n";

        csv << "longline_decomposition," << obs.year << "," << (a + 1) << ","
            << n[i] << "," << sel_longline[i] << "," << selected[i] << ","
            << pred << "," << obs_a << "," << resid << "," << abs_resid << ","
            << pred << "," << obs.fleet << "," << obs.index << ","
            << obs.catch_mt << "\n";

        residual_rows.push_back({obs.year, a + 1, obs_a, pred, resid, abs_resid,
                                 n[i], sel_longline[i], selected[i]});
      }
    }

    std::array<double, kAges> total_sel{};
    for (int a = 0; a < kAges; ++a) {
      const std::size_t i = static_cast<std::size_t>(a);
      total_sel[i] = sel_longline[i] + sel_purse_seine[i];
    }

    std::array<double, kAges> next{};
    const auto weight = default_weight_at_age();
    const auto maturity = default_maturity_at_age();
    const double spawning_biomass =
        level23_bh_sync::spawning_biomass_from_numbers(n, weight, maturity);
    const double expected_recruitment =
        level23_bh_sync::beverton_holt_recruitment(spawning_biomass, r0,
                                                   steepness, phi0);
    next[0] = expected_recruitment *
              std::exp(fit.u_hat[static_cast<Eigen::Index>(t)]);
    for (int a = 1; a < kAges; ++a) {
      const std::size_t prev = static_cast<std::size_t>(a - 1);
      const double z_prev = m + fbar * total_sel[prev];
      next[static_cast<std::size_t>(a)] = n[prev] * std::exp(-z_prev);
    }
    const std::size_t last = static_cast<std::size_t>(kAges - 1);
    const double z_last = m + fbar * total_sel[last];
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
         "longline_selectivity,selected_numbers\n";

  for (std::size_t i = 0; i < std::min<std::size_t>(25, residual_rows.size());
       ++i) {
    const auto &r = residual_rows[i];
    txt << r.year << "," << r.age << "," << r.observed << "," << r.predicted
        << "," << r.residual << "," << r.abs_residual << "," << r.n_at_age
        << "," << r.sel << "," << r.selected << "\n";
  }

  std::array<double, kAges> mean_obs{};
  std::array<double, kAges> mean_pred{};
  std::array<double, kAges> mean_resid{};
  std::array<double, kAges> mean_abs{};
  std::array<int, kAges> counts{};

  for (const auto &r : residual_rows) {
    const std::size_t i = static_cast<std::size_t>(r.age - 1);
    mean_obs[i] += r.observed;
    mean_pred[i] += r.predicted;
    mean_resid[i] += r.residual;
    mean_abs[i] += r.abs_residual;
    counts[i] += 1;
  }

  txt << "\nMean by age\n";
  txt << "-----------\n";
  txt << "age,mean_observed,mean_predicted,mean_residual,mean_abs_residual,n\n";
  for (int a = 0; a < kAges; ++a) {
    const std::size_t i = static_cast<std::size_t>(a);
    if (counts[i] > 0) {
      mean_obs[i] /= static_cast<double>(counts[i]);
      mean_pred[i] /= static_cast<double>(counts[i]);
      mean_resid[i] /= static_cast<double>(counts[i]);
      mean_abs[i] /= static_cast<double>(counts[i]);
    }
    txt << (a + 1) << "," << mean_obs[i] << "," << mean_pred[i] << ","
        << mean_resid[i] << "," << mean_abs[i] << "," << counts[i] << "\n";
    csv << "mean_by_age,," << (a + 1) << ",,,," << mean_pred[i] << ","
        << mean_obs[i] << "," << mean_resid[i] << "," << mean_abs[i]
        << ",,longline,,\n";
  }
}

} // namespace pifsc_bigeye_tuna

using pifsc_bigeye_tuna::write_level18_longline_prediction_decomposition;
