/**
 * @file report.hpp
 * @brief Historical reporting with explicit units and conventions.
 *
 * Reports use numbers in thousands, biomass/catch in tonnes and pooled SS3
 * units-4 F. The caller supplies an existing output directory. Files here are
 * scratch artifacts; the supported release runner isolates them and publishes
 * the entire completed bundle atomically after validation, rather than treating
 * one CSV as a release.
 */
#pragma once
#include "benchmarks.hpp"
#include "likelihood.hpp"
#include <filesystem>
#include <fstream>
#include <iomanip>

namespace swordfish {
/// Write independently calculated reference points for this fitted model state.
/// The standard table includes unfished, Btgt20, SPR20 and MSY. Biomass/catch
/// are tonnes, recruitment is thousands, and annual F uses SS3 units 4 ages
/// 1..10. The output directory must already exist. Failures propagate to the
/// fitting caller; SS3 values are used only by separate comparison checks.
inline void write_reference_points(const InputData &data,
                                   const ParameterValues<double> &p,
                                   const Deviations<double> &dev,
                                   const std::filesystem::path &folder) {
  const auto points = management_reference_points(data, p, dev);
  std::ofstream out(folder / "management_reference_points.csv");
  out << "benchmark,F_multiplier,SPR,SSB_tonnes,recruitment_thousands,dead_"
         "catch_tonnes,annual_F_SS3_units4_ages1_10\n"
      << std::setprecision(16);
  for (const auto &point : points) {
    const auto &r = point.second;
    if (!std::isfinite(r.scale) || !std::isfinite(r.spr) ||
        !std::isfinite(r.ssb) || !std::isfinite(r.recruitment) ||
        !std::isfinite(r.yield) || !std::isfinite(r.annual_f))
      throw std::runtime_error("Nonfinite model reference point");
    out << point.first << ',' << r.scale << ',' << r.spr << ',' << r.ssb << ','
        << r.recruitment << ',' << r.yield << ',' << r.annual_f << '\n';
  }
  out.flush();
  if (!out)
    throw std::runtime_error("Failed writing model reference points");
}
/// Write historical annual/seasonal trajectories, fitted indices, length
/// diagnostics, likelihood components and management reference points to the
/// supplied output folder. Release publication is performed by the isolated
/// runner.
/// @param data Validated original-input schema for the supported swordfish
/// assessment.
/// @param p Native canonical-name parameter map; all required fixed and
/// estimated values must be present.
/// @param dev Year-to-log-recruitment-deviation map; the caller enforces the
/// main-period constraint.
/// @param folder Existing scratch output directory; publication is managed by
/// the release runner.
inline void write_fitted_outputs(const InputData &data,
                                 const ParameterValues<double> &p,
                                 const Deviations<double> &dev,
                                 const std::filesystem::path &folder) {
  const auto run = run_population(data, p, dev);
  const auto fit = evaluate_likelihood(data, p, dev, run);
  std::ofstream seasonal(folder / "quadra_fitted_seasonal.csv"),
      annual(folder / "quadra_fitted_annual.csv"),
      indices(folder / "quadra_fitted_index_fits.csv"),
      lengths(folder / "quadra_fitted_length_fits.csv"),
      components(folder / "quadra_fitted_likelihood.csv");
  if (!seasonal || !annual || !indices || !lengths || !components)
    throw std::runtime_error("Cannot write fitted assessment outputs");
  seasonal << "year,season,biomass,age1plus_biomass,spawning_biomass,"
              "recruitment_at_settlement";
  for (int f = 1; f <= nfishing; ++f)
    seasonal << ",apical_F_" << f;
  seasonal << '\n' << std::setprecision(15);
  for (const auto &r : run.seasons) {
    seasonal << r.year << ',' << r.season << ',' << r.biomass << ','
             << r.summary_biomass << ',';
    if (r.ssb >= 0)
      seasonal << r.ssb;
    seasonal << ',' << (r.season == 3 ? r.recruits : 0);
    for (double f : r.hybrid.apical_f)
      seasonal << ',' << f;
    seasonal << '\n';
  }
  annual << "year,biomass_january,age1plus_biomass_january,spawning_biomass,"
            "recruitment_thousands,annual_F_SS3_units4_ages1_10,dead_catch_"
            "tonnes\n"
         << std::setprecision(15);
  for (int year = data.dimensions.first_year; year <= data.dimensions.last_year;
       ++year) {
    const auto &start = population_row(run, data, year, 1);
    const auto &end = population_row(run, data, year, 4);
    double natural_end = 0, actual_end = 0, catch_biomass = 0;
    for (int sex = 0; sex < 2; ++sex)
      for (int age = 1; age <= 10; ++age) {
        const int j = cell(sex, age);
        natural_end +=
            start.beginning[j] * std::exp(-data.biology.natural_mortality[j]);
        actual_end += end.removals.survivors[j];
      }
    const auto curves = length_selectivity(data, p, year);
    for (int s = 1; s <= 4; ++s) {
      const auto &row = population_row(run, data, year, s);
      const auto biology = seasonal_biology(data, p, curves, s);
      for (int f = 0; f < nfishing; ++f)
        for (int j = 0; j < 32; ++j)
          catch_biomass += row.removals.catch_numbers[f][j] *
                           biology.selected_biomass[f][j] /
                           biology.selectivity[f][j];
    }
    annual << year << ',' << start.biomass << ',' << start.summary_biomass
           << ',' << population_row(run, data, year, 2).ssb << ','
           << population_row(run, data, year, 3).recruits << ','
           << std::log(natural_end / actual_end) << ',' << catch_biomass
           << '\n';
  }
  indices << "year,season,fleet,observed,expected,log_q,se,log_residual,nll,"
             "lambda\n"
          << std::setprecision(15);
  for (const auto &r : fit.indices)
    indices << r.year << ',' << r.season << ',' << r.fleet << ','
            << r.observation << ',' << r.prediction << ',' << r.log_q << ','
            << r.se << ',' << r.deviation << ',' << r.nll << ','
            << likelihood_weight(data, 1, r.fleet, 5) << '\n';
  lengths << "year,season,fleet,sample_size,observed_mean,expected_mean,nll,"
             "lambda\n"
          << std::setprecision(15);
  for (const auto &r : fit.lengths)
    lengths << r.year << ',' << r.season << ',' << r.fleet << ','
            << r.sample_size << ',' << r.observed_mean << ',' << r.expected_mean
            << ',' << r.nll << ',' << likelihood_weight(data, 4, r.fleet, 5)
            << '\n';
  components << "component,value\n" << std::setprecision(15);
  for (const auto &r : fit.components)
    components << r.first << ',' << r.second << '\n';
  components << "TOTAL," << fit.total << '\n';
  if (!seasonal || !annual || !indices || !lengths || !components)
    throw std::runtime_error("Failed writing fitted assessment outputs");
  write_reference_points(data, p, dev, folder);
}
} // namespace swordfish
