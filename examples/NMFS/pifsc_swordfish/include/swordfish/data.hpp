/**
 * @file data.hpp
 * @brief Population schema and canonical index conventions.
 *
 * State cells are sex-major: female=0, male=1, ages 0..15 with age 15
 * accumulating older fish. Fleet identifiers in source observations are
 * one-based. Counts use thousands of fish; individual weights use kg, so
 * biomass and dead yield are tonnes. These compile-time dimensions describe the
 * supplied assessment, not arbitrary SS3 models.
 */
#pragma once
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

namespace swordfish {
constexpr int nsex = 2;
constexpr int nage = 16;
constexpr int nfishing = 19;
/// Fixed 2-sex by 16-age abundance array. Cells are thousands of fish; age 15
/// is the plus group.
/// @details Index with cell(sex,age), not a one-based SS3 sex/age code. Use
/// zero_state before accumulation: default construction of supported AD types
/// is not a numerical-zero contract.
template <class T> using State = std::array<T, nsex * nage>;
// Quadra AD scalars have an intentionally uninitialized default constructor.
// Value-initializing an array is insufficient; construct numerical zeros
// explicitly.
/// Construct an explicitly zero-filled sex-by-age state, including AD scalars
/// whose default constructor does not initialize a numerical zero.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active
/// evaluation tape.
/// @return A new state containing numerical T(0) in every cell.
template <class T> State<T> zero_state() {
  State<T> state;
  state.fill(T(0));
  return state;
}
/// Map zero-based sex and age to the flattened 32-cell population-state index;
/// callers supply sex 0-1 and age 0-15.
/// @param sex Zero-based sex index: 0 female, 1 male.
/// @param age Integer model age 0..15, with age 15 representing the plus group.
/// @return Sex-major offset sex*nage+age; no bounds check is performed by this
/// primitive.
constexpr int cell(int sex, int age) { return sex * nage + age; }
/// Assessment calendar, four seasonal durations, spawning month, area count and
/// plus-group age.
/// @details first_year/last_year are inclusive calendar years.
/// seasons/subseasons and plus_age describe the supported discretization;
/// season_months contains durations, not start months. spawn_month is one-based
/// calendar timing.
struct Dimensions {
  int first_year, last_year, seasons, subseasons, plus_age, areas;
  double spawn_month;
  std::array<double, 4> season_months;
};
/// One-based fleet identity, type, observation timing, catch units, multiplier
/// flag and source name.
/// @details id is one-based and matches observation records. type distinguishes
/// fishing/index fleets; catch_units=1 is biomass and 2 is numbers. timing
/// controls point versus integrated sampling. This schema has one area and no
/// supported catch multipliers.
struct Fleet {
  int id, type, timing, area, catch_units, catch_multiplier;
  std::string name;
};
/// One observed fleet-season catch with year, value and log-scale standard
/// error; fleet units select tonnes or thousands of fish.
/// @details value follows the referenced fleet units, while log_se is a
/// log-scale error. year=-999 carries the supported initialization observation;
/// normal rows are historical year/season catches, not estimates of apical F.
struct Catch {
  int year, season, fleet;
  double value, log_se;
};
/// Abundance-index observation and inclusion flag with source month, fleet,
/// positive observation and log standard error.
/// @details month retains source fractional timing; value and log_se must be
/// positive for included lognormal observations. use controls inclusion without
/// deleting provenance rows. Special recruitment-unit indices are handled
/// explicitly in the likelihood.
struct Index {
  int year;
  double month;
  int fleet;
  double value, log_se;
  bool use;
};
/// Length-composition sample, inclusion flag, sex/partition coding, sample size
/// and both source blocks of 51 bins.
/// @details values retains both 51-bin source blocks. For the supplied
/// combined-sex data the second block duplicates the first and is not a second
/// independent sample. Negative source years set use=false; sample_size already
/// contains the supplied weighting.
struct Composition {
  int year;
  double month;
  int fleet, sex, partition;
  double sample_size;
  bool use;
  std::array<double, 102> values;
};
/// Named SS3 control specification, original numeric columns and source line
/// for provenance.
/// @details columns preserves the source control order: lower, upper, initial
/// and estimation metadata. Estimated phases are read from column 6; extended
/// selectivity rows also carry block number/function in columns 12/13.
/// source_line supports diagnostics.
struct Parameter {
  std::string name;
  std::vector<double> columns;
  int source_line;
};
/// Sex-specific von Bertalanffy reference ages/lengths, coefficient of
/// variation and length-weight coefficients.
/// @details Lengths are cm and k is annual. Reference lengths are tied to
/// age1/age2; cv_young/cv_old are coefficients of variation, not length SDs.
/// weight_coefficient and weight_power map individual cm length to kg weight.
struct Growth {
  double age1, age2, length1, length2, k, cv_young, cv_old;
  double weight_coefficient, weight_power;
};
/// Fixed natural mortality, sex-specific growth, maturity, sex allocation and
/// settlement settings.
/// @details natural_mortality is sex-major annual M. growth has one record per
/// sex; maturity and settlement settings describe the supplied control setup.
/// Fitted native parameter maps remain authoritative for parameters evaluated
/// by the biological functions.
struct Biology {
  State<double> natural_mortality;
  std::array<Growth, 2> growth;
  double maturity_l50, maturity_slope, female_fraction;
  int first_mature_age;
  double settlement_month, settlement_age;
};
/// Complete supported swordfish observation and control schema; unsupported SS3
/// variants are rejected rather than approximated silently.
/// @details This object owns parsed observations and controls; no trajectory
/// from Report.sso is a population driver. Arrays preserve source codes and
/// ordering, and validate checks the explicitly supported schema. Further
/// method-specific layout checks occur before likelihood/projection evaluation.
struct InputData {
  Dimensions dimensions;
  Biology biology;
  std::vector<Fleet> fleets;
  std::vector<Catch> catches;
  std::vector<Index> indices;
  std::vector<Composition> compositions;
  std::vector<double> population_length_bins, length_bins;
  std::vector<std::array<int, 4>> index_specs, size_patterns, age_patterns;
  std::vector<std::array<double, 7>> length_controls;
  std::vector<Parameter> parameters;
  std::vector<std::vector<double>> blocks;
  std::vector<std::array<double, 6>> q_setup;
  std::vector<std::array<double, 3>> variance_adjustments;
  std::vector<std::array<double, 5>> lambda_changes;
  int recdev_start, recdev_end, recdev_phase, f_method, hybrid_iterations,
      max_lambda_phase, sd_offset;
  double max_f, plus_group_growth;
  std::array<double, 13> recruitment_advanced;

  /// Check supported dimensions, season durations, mortality, fleet/catch
  /// definitions and growth before evaluating the imported data.
  void validate() const {
    if (dimensions.first_year > dimensions.last_year ||
        dimensions.seasons != 4 || dimensions.plus_age != 15 ||
        dimensions.areas != 1 || fleets.size() != 27)
      throw std::invalid_argument("Unsupported swordfish dimensions");
    double months = 0;
    for (double duration : dimensions.season_months) {
      if (!(duration > 0) || !std::isfinite(duration))
        throw std::invalid_argument("Season duration");
      months += duration;
    }
    if (std::abs(months - 12) > 1e-10)
      throw std::invalid_argument("Seasons must span one year");
    for (double m : biology.natural_mortality)
      if (!(m >= 0) || !std::isfinite(m))
        throw std::invalid_argument("Natural mortality");
    if (!std::isfinite(biology.female_fraction) ||
        biology.female_fraction < 0 || biology.female_fraction > 1)
      throw std::invalid_argument("Sex allocation");
    for (const auto &f : fleets)
      if (f.id < 1 || f.id > 27 || (f.type != 1 && f.type != 3) ||
          (f.catch_units != 1 && f.catch_units != 2))
        throw std::invalid_argument("Fleet definition");
    for (const auto &c : catches)
      if (c.fleet < 1 || c.fleet > nfishing || c.season < 1 || c.season > 4 ||
          !std::isfinite(c.value) || c.value < 0 || !std::isfinite(c.log_se) ||
          c.log_se < 0)
        throw std::invalid_argument("Catch observation");
    if (population_length_bins.empty())
      throw std::invalid_argument("Missing length bins");
    for (const auto &g : biology.growth)
      if (!(g.age2 > g.age1 && g.age1 > 0 && g.k > 0 && g.length2 > g.length1 &&
            g.length1 > population_length_bins.front()))
        throw std::invalid_argument("Growth definition");
  }
};
} // namespace swordfish
