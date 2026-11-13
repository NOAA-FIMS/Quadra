/**
 * @file dynamics.hpp
 * @brief Seasonal mortality, calendar aging and recruitment settlement.
 *
 * Mortality inputs are annual instantaneous rates and exposure times are years. All abundance and
 * removals share the thousands-of-fish scale. Callers provide nonnegative abundance/rates. Sex
 * allocation and January plus-group aging remain separate from mortality so July settlement is not
 * inadvertently aged when a season advances.
 */
#pragma once
#include "data.hpp"
#include <cmath>
#include <vector>

namespace swordfish {
// All rates are annual instantaneous rates; duration and elapsed are in years.
// Counts are thousands of fish. Multiplying by kg/fish yields metric tonnes.
/// Calculate the stable Baranov exposure (1-exp(-Z*duration))/Z, using a near-zero expansion; Z is an annual rate and duration is years.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active evaluation tape.
/// @param z Total annual instantaneous mortality; callers supply nonnegative rates.
/// @param duration Exposure duration in calendar years; a quarter is 0.25.
/// @return Integrated survival per unit abundance; its limit as z tends to zero is duration.
template <class T> T exposure(const T& z, double duration) {
  using std::exp;
  const T x = z * duration;
  // Finite derivative and exact limit at zero mortality, including under AD.
  if (x < 1e-5)
    return duration * (T(1) - x/T(2) + x*x/T(6) - x*x*x/T(24) + x*x*x*x/T(120));
  return (T(1) - exp(-x)) / z;
}
/// Seasonal survivors, natural deaths and fleet catch numbers, all on the same thousands-of-fish scale as the starting state.
/// @details survivors and natural_deaths are sex-major states; catch_numbers is indexed by fleet then sex/age. advance_season writes every returned cell. Numbers are thousands, so population mass balance is tested on that common scale.
template <class T> struct SeasonResult {
  State<T> survivors{}, natural_deaths{};
  std::vector<State<T>> catch_numbers;
};
/// Apply natural and fleet fishing mortality over a fraction of a year, returning survivors and separately attributed removals without altering the input state.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active evaluation tape.
/// @param numbers Sex-by-age abundance at the start of the calculation, in thousands of fish.
/// @param m Fixed sex-by-age natural mortality, expressed as annual instantaneous rates.
/// @param fleet_f Per-fleet sex-by-age fishing mortality, expressed as annual instantaneous rates.
/// @param duration Exposure duration in calendar years; a quarter is 0.25.
/// @return Survivors, natural deaths and fleet catch numbers on the input abundance scale.
template <class T> SeasonResult<T> advance_season(
    const State<T>& numbers, const State<double>& m,
    const std::vector<State<T>>& fleet_f, double duration) {
  using std::exp;
  if (!(duration > 0 && duration <= 1)) throw std::invalid_argument("Season duration");
  SeasonResult<T> result;
  result.catch_numbers.resize(fleet_f.size());
  for (int j = 0; j < nsex*nage; ++j) {
    T z = T(m[j]);
    for (const auto& f : fleet_f) z += f[j];
    const T exposed = numbers[j] * exposure(z, duration);
    result.survivors[j] = numbers[j] * exp(-z * duration);
    result.natural_deaths[j] = exposed * m[j];
    for (size_t f = 0; f < fleet_f.size(); ++f)
      result.catch_numbers[f][j] = exposed * fleet_f[f][j];
  }
  return result;
}
/// Advance survivors one age at the calendar-year boundary, retain accumulated ages in the plus group and clear age zero.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active evaluation tape.
/// @param survivors Post-mortality abundance in thousands of fish, before January aging.
/// @return A new January state, retaining old plus-group survivors and admitting age-14 survivors.
template <class T> State<T> age_on_january_first(const State<T>& survivors) {
  State<T> next=zero_state<T>();
  for (int s = 0; s < nsex; ++s) {
    for (int a = 1; a < nage; ++a) next[cell(s,a)] = survivors[cell(s,a-1)];
    next[cell(s,nage-1)] += survivors[cell(s,nage-1)];
  }
  return next;
}
/// Mutate age-zero abundance by adding recruitment allocated to females and males using the supplied female fraction.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active evaluation tape.
/// @param numbers Sex-by-age abundance at the start of the calculation, in thousands of fish.
/// @param recruits Recruitment in thousands of fish to allocate at the July settlement event.
/// @param female_fraction Recruitment fraction assigned to females; male allocation is its complement.
/// @post Only age-zero cells are incremented; other ages and caller parameter values are unchanged.
template <class T, class Fraction> void settle(State<T>& numbers, const T& recruits,
                               Fraction female_fraction) {
  numbers[cell(0,0)] += recruits * female_fraction;
  numbers[cell(1,0)] += recruits * (1-female_fraction);
}
/// Calculate female mature biomass after survival to a calendar-correct within-season spawning time.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active evaluation tape.
/// @param numbers Sex-by-age abundance at the start of the calculation, in thousands of fish.
/// @param m Fixed sex-by-age natural mortality, expressed as annual instantaneous rates.
/// @param fleet_f Per-fleet sex-by-age fishing mortality, expressed as annual instantaneous rates.
/// @param mature_weight Mature female weight contribution per fish in kg, indexed by age.
/// @param elapsed Calendar-year fraction between seasonal beginning and spawning.
/// @return Female spawning biomass in tonnes under calendar-correct survival.
template <class T> T spawning_biomass(const State<T>& numbers,
    const State<double>& m, const std::vector<State<T>>& fleet_f,
    const std::array<double,nage>& mature_weight, double elapsed) {
  using std::exp;
  if (!(elapsed >= 0 && elapsed <= 1)) throw std::invalid_argument("Spawn timing");
  T ssb = T(0);
  for (int a = 0; a < nage; ++a) {
    T z = T(m[cell(0,a)]);
    for (const auto& f : fleet_f) z += f[cell(0,a)];
    ssb += numbers[cell(0,a)] * exp(-z*elapsed) * mature_weight[a];
  }
  return ssb;
}
// Version-specific reproduction: SS3 3.30.20 multiplies the annual Z rate
// directly by the fraction of the spawning season elapsed (SS_popdyn.tpl).
// Keep this separate from the calendar-year convention above. This is a
// compatibility calculation, not an assertion that the timings are equivalent.
/// Reproduce the legacy SS3 v3.30.20 spawning-survival timing convention used by the supplied benchmark, deliberately distinct from the calendar-correct helper.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active evaluation tape.
/// @param numbers Sex-by-age abundance at the start of the calculation, in thousands of fish.
/// @param m Fixed sex-by-age natural mortality, expressed as annual instantaneous rates.
/// @param fleet_f Per-fleet sex-by-age fishing mortality, expressed as annual instantaneous rates.
/// @param mature_weight Mature female weight contribution per fish in kg, indexed by age.
/// @param calendar_elapsed Calendar-year fraction from seasonal beginning to spawning.
/// @param season_duration Season length in years; used to reconstruct the legacy SS3 within-season fraction.
/// @return Female spawning biomass in tonnes under legacy v3.30.20 survival.
/// @warning calendar_elapsed/season_duration is used directly against annual mortality to reproduce this SS3 version. For May spawning in Q2, this is 1/3 rather than 1/12 calendar years; use the separate calendar helper for that alternative definition.
template <class T> T spawning_biomass_ss33020(const State<T>& numbers,
    const State<double>& m, const std::vector<State<T>>& fleet_f,
    const std::array<double,nage>& mature_weight,
    double calendar_elapsed, double season_duration) {
  if (!(season_duration > 0 && calendar_elapsed >= 0 && calendar_elapsed <= season_duration))
    throw std::invalid_argument("SS3 spawning timing");
  return spawning_biomass(numbers,m,fleet_f,mature_weight,calendar_elapsed/season_duration);
}
// Real age is elapsed time since birth. The aggregate plus group's mean size
// needs SS3's additional growth adjustment and is not supplied by this function.
/// Evaluate the ordinary von Bertalanffy growth curve at a real age using fixed growth settings and length at birth.
/// @param real_age Elapsed years since birth, including fractional age.
/// @param g Fixed growth reference ages, lengths, rate and variability for one sex.
/// @param length_at_birth Mean length at birth in cm.
/// @return The typed result described above, with the module units and indexing conventions.
inline double mean_length(double real_age, const Growth& g, double length_at_birth) {
  if (real_age <= 0) return length_at_birth;
  if (real_age < g.age1)
    return length_at_birth + (g.length1-length_at_birth)*real_age/g.age1;
  return g.length1 + (g.length2-g.length1) *
      (-std::expm1(-g.k*(real_age-g.age1))) / (-std::expm1(-g.k*(g.age2-g.age1)));
}
/// Return expected recruitment from spawning biomass, virgin recruitment, virgin spawning biomass and steepness; retain scalar-generic derivatives.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active evaluation tape.
/// @param ssb Female spawning biomass in tonnes under the specified timing convention.
/// @param r0 Virgin expected annual recruitment in thousands of fish.
/// @param b0 Virgin spawning biomass in tonnes, consistent with r0.
/// @param h Beverton-Holt steepness on its native scale.
/// @return Expected recruitment in thousands of fish; no stochastic deviation is drawn.
template <class T, class H> T beverton_holt(const T& ssb, const T& r0, const T& b0, H h) {
  return T(4*h) * r0 * ssb / (b0*T(1-h) + ssb*T(5*h-1));
}
} // namespace swordfish
