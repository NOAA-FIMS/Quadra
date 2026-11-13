/**
 * @file biology.hpp
 * @brief Differentiable growth, length integration and selectivity.
 *
 * Lengths are cm, weights kg per fish, sex is zero-based, and within-year times
 * are fractions of a calendar year. Parameters are native canonical-name
 * values, not optimizer coordinates. T may be double or a Quadra AD scalar;
 * cache only within the active evaluation/tape. Unsupported selectivity
 * patterns and block functions must fail rather than silently approximate
 * another model.
 */
#pragma once
#include "dynamics.hpp"
#include <algorithm>
#include <map>
#include <string>

namespace swordfish {
/// Map canonical parameter names to native values of a common scalar type.
/// @details Keys are canonical control names. Values are native model values,
/// not logits; all entries used together must share the same scalar
/// type/evaluation tape. Missing required keys fail through checked map lookup.
template <class T> using ParameterValues = std::map<std::string, T>;
/// Populate a parameter map from control initial values; these may differ from
/// fitted Report.sso values.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active
/// evaluation tape.
/// @param data Parsed controls whose initial native values are read;
/// report-fitted values are separate fixtures.
/// @return The typed result described above, with the module units and indexing
/// conventions.
template <class T> ParameterValues<T> initial_values(const InputData &data) {
  ParameterValues<T> p;
  for (const auto &spec : data.parameters)
    p.emplace(spec.name, T(spec.columns[2]));
  return p;
}
/// Calculate a numerically stable inverse-logit for ordinary or AD scalars.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active
/// evaluation tape.
/// @param x Unconstrained scalar log-odds; evaluated with a value-based stable
/// branch.
/// @return Inverse-logit probability on the same scalar/AD tape as x.
template <class T> T logistic(const T &x) {
  using std::exp;
  if (x < 0) {
    const T e = exp(x);
    return e / (T(1) + e);
  }
  return T(1) / (T(1) + exp(-x));
}
// Standard normal CDF approximation (maximum absolute error < 7.5e-8).
// Written in scalar arithmetic so the length integration remains
// differentiable.
/// Approximate the standard normal CDF with scalar-generic operations, enabling
/// differentiation through age-length probabilities.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active
/// evaluation tape.
/// @param x Standard-normal deviate; derivatives follow the smooth
/// approximation.
/// @return Approximate normal cumulative probability on the same scalar/AD tape
/// as x.
template <class T> T normal_cdf(const T &x) {
  using std::exp;
  const bool negative = x < 0;
  const T a = negative ? -x : x;
  const T t = T(1) / (T(1) + T(0.2316419) * a);
  const T tail =
      T(0.3989422804014327) * exp(-a * a / T(2)) * t *
      (T(0.319381530) +
       t * (T(-0.356563782) +
            t * (T(1.781477937) + t * (T(-1.821255978) + t * T(1.330274429)))));
  return negative ? tail : T(1) - tail;
}
/// Integrate normal lengths across lower-edge bins, accumulating the two tails
/// in the end bins and rejecting invalid grids or scale.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active
/// evaluation tape.
/// @param edges Strictly increasing lower length-bin edges in cm; end bins
/// accumulate the tails.
/// @param mean Mean individual length in cm.
/// @param sd Strictly positive length standard deviation in cm.
/// @return The typed result described above, with the module units and indexing
/// conventions.
template <class T>
std::vector<T> length_probabilities(const std::vector<double> &edges,
                                    const T &mean, const T &sd) {
  if (edges.size() < 2 || !(sd > 0))
    throw std::invalid_argument("Length distribution inputs");
  std::vector<T> probability(edges.size());
  T previous = T(0);
  for (size_t b = 0; b + 1 < edges.size(); ++b) {
    const T cumulative = normal_cdf((T(edges[b + 1]) - mean) / sd);
    probability[b] = cumulative - previous;
    previous = cumulative;
  }
  probability.back() =
      T(1) - previous; // both open tails absorbed into end bins
  return probability;
}
/// Calculate sex-specific mean length with the supplied birth timing, seasonal
/// growth and SS3 plus-group growth adjustment.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active
/// evaluation tape.
/// @param data Validated original-input schema for the supported swordfish
/// assessment.
/// @param p Native canonical-name parameter map; all required fixed and
/// estimated values must be present.
/// @param sex Zero-based sex index: 0 female, 1 male.
/// @param age Integer model age 0..15, with age 15 representing the plus group.
/// @param year_fraction Time within the calendar year in years; not the
/// fraction of the current season.
/// @return Mean length in cm, including the supported plus-group adjustment.
template <class T>
T length_at_age(const InputData &data, const ParameterValues<T> &p, int sex,
                int age, double year_fraction) {
  using std::exp;
  if (sex < 0 || sex >= nsex || age < 0 || age >= nage || year_fraction < 0 ||
      year_fraction >= 1)
    throw std::invalid_argument("Biological age/time index");
  const std::string suffix = sex == 0 ? "_Fem_GP_1" : "_Mal_GP_1";
  const T l1 = p.at("L_at_Amin" + suffix), l2 = p.at("L_at_Amax" + suffix),
          k = p.at("VonBert_K" + suffix);
  const auto &g = data.biology.growth[sex];
  const T e = exp(-k * (g.age2 - g.age1));
  const T linf = (l2 - l1 * e) / (T(1) - e);
  const double birth = (data.biology.settlement_month - 1) / 12;
  auto ordinary = [&](double a) -> T {
    if (a <= 0)
      return T(data.population_length_bins.front());
    if (a < g.age1)
      return T(data.population_length_bins.front()) +
             (l1 - T(data.population_length_bins.front())) * (a / g.age1);
    return linf + (l1 - linf) * exp(-k * (a - g.age1));
  };
  if (age < nage - 1)
    return ordinary(age + year_fraction - birth);
  if (data.plus_group_growth != -999)
    throw std::invalid_argument("Unsupported plus-group growth option");
  const T exact = ordinary(age - birth);
  T numerator = T(0);
  double denominator = 0;
  for (int j = 0; j <= data.dimensions.plus_age; ++j) {
    const double w = std::exp(-0.2 * j);
    numerator += T(w) * (exact + (linf - exact) *
                                     (double(j) / data.dimensions.plus_age));
    denominator += w;
  }
  const T january = numerator / denominator;
  return linf + (january - linf) * exp(-k * year_fraction);
}
/// Interpolate age-dependent length variability and convert its coefficient of
/// variation to a standard deviation.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active
/// evaluation tape.
/// @param data Validated original-input schema for the supported swordfish
/// assessment.
/// @param p Native canonical-name parameter map; all required fixed and
/// estimated values must be present.
/// @param sex Zero-based sex index: 0 female, 1 male.
/// @param age Integer model age 0..15, with age 15 representing the plus group.
/// @param time Within-year time passed through to the mean-length calculation.
/// @param mean Mean individual length in cm.
/// @return Length standard deviation in cm.
template <class T>
T length_sd(const InputData &data, const ParameterValues<T> &p, int sex,
            int age, double time, const T &mean) {
  const std::string suffix = sex == 0 ? "_Fem_GP_1" : "_Mal_GP_1";
  const T young = p.at("CV_young" + suffix), old = p.at("CV_old" + suffix);
  const double real_age = age + time - (data.biology.settlement_month - 1) / 12;
  const auto &g = data.biology.growth[sex];
  T cv = young;
  if (real_age >= g.age2)
    cv = old;
  else if (real_age >= g.age1)
    cv = young + (old - young) * (mean - p.at("L_at_Amin" + suffix)) /
                     (p.at("L_at_Amax" + suffix) - p.at("L_at_Amin" + suffix));
  return mean * cv;
}
/// Return bin probabilities for one sex, age and within-year time using fitted
/// growth parameters.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active
/// evaluation tape.
/// @param data Validated original-input schema for the supported swordfish
/// assessment.
/// @param p Native canonical-name parameter map; all required fixed and
/// estimated values must be present.
/// @param sex Zero-based sex index: 0 female, 1 male.
/// @param age Integer model age 0..15, with age 15 representing the plus group.
/// @param time Within-year time in years at which growth/length integration is
/// evaluated.
/// @return The typed result described above, with the module units and indexing
/// conventions.
template <class T>
std::vector<T> age_length_key(const InputData &data,
                              const ParameterValues<T> &p, int sex, int age,
                              double time) {
  const T mean = length_at_age(data, p, sex, age, time);
  return length_probabilities(data.population_length_bins, mean,
                              length_sd(data, p, sex, age, time, mean));
}
/// Map a population length-bin lower edge to its integration midpoint,
/// including the final bin convention.
/// @param data Validated original-input schema for the supported swordfish
/// assessment.
/// @param bin Zero-based index into the population length integration grid.
/// @return Integration midpoint in cm.
inline double length_midpoint(const InputData &data, size_t bin) {
  const auto &e = data.population_length_bins;
  const double width =
      bin + 1 < e.size() ? e[bin + 1] - e[bin] : e.back() - e[e.size() - 2];
  return e[bin] + width / 2;
}
/// Evaluate the sex-specific length-weight power relationship; length is cm and
/// individual weight is kg.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active
/// evaluation tape.
/// @param p Native canonical-name parameter map; all required fixed and
/// estimated values must be present.
/// @param sex Zero-based sex index: 0 female, 1 male.
/// @param length Individual length in cm.
/// @return Individual weight in kg.
template <class T>
T weight_at_length(const ParameterValues<T> &p, int sex, double length) {
  using std::exp;
  const std::string suffix = sex == 0 ? "_Fem_GP_1" : "_Mal_GP_1";
  return p.at("Wtlen_1" + suffix) *
         exp(p.at("Wtlen_2" + suffix) * std::log(length));
}
/// Evaluate the length-based logistic female maturity curve from its midpoint
/// and slope.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active
/// evaluation tape.
/// @param p Native canonical-name parameter map; all required fixed and
/// estimated values must be present.
/// @param length Individual length in cm.
/// @return Mature fraction at the supplied length.
template <class T>
T maturity_at_length(const ParameterValues<T> &p, double length) {
  return logistic(-p.at("Mat_slope_Fem_GP_1") *
                  (T(length) - p.at("Mat50%_Fem_GP_1")));
}
/// Integrate individual weight over the age-length key at a supplied
/// within-year time.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active
/// evaluation tape.
/// @param data Validated original-input schema for the supported swordfish
/// assessment.
/// @param p Native canonical-name parameter map; all required fixed and
/// estimated values must be present.
/// @param sex Zero-based sex index: 0 female, 1 male.
/// @param age Integer model age 0..15, with age 15 representing the plus group.
/// @param time Within-year time in years at which growth/length integration is
/// evaluated.
/// @return Mean individual weight in kg.
template <class T>
T weight_at_age(const InputData &data, const ParameterValues<T> &p, int sex,
                int age, double time) {
  const auto key = age_length_key(data, p, sex, age, time);
  T weight = T(0);
  for (size_t b = 0; b < key.size(); ++b)
    weight += key[b] * weight_at_length(p, sex, length_midpoint(data, b));
  return weight;
}
// SS3 v3.30.20 uses the start of the subseason containing the spawn month
// for the maturity ALK. With these inputs that is the beginning of quarter 2.
/// Locate the start of the subseason containing spawning, matching the supplied
/// SS3 maturity integration convention.
/// @param data Validated original-input schema for the supported swordfish
/// assessment.
/// @return The typed result described above, with the module units and indexing
/// conventions.
inline double maturity_integration_time(const InputData &data) {
  const double spawn = (data.dimensions.spawn_month - 1) / 12;
  double start = 0;
  for (double months : data.dimensions.season_months) {
    const double duration = months / 12;
    if (spawn < start + duration) {
      const double step = duration / data.dimensions.subseasons;
      return start + std::floor((spawn - start) / step + 1e-9) * step;
    }
    start += duration;
  }
  throw std::invalid_argument("Spawn month outside year");
}
/// Integrate mature female biomass contribution across length bins, respecting
/// the first mature age.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active
/// evaluation tape.
/// @param data Validated original-input schema for the supported swordfish
/// assessment.
/// @param p Native canonical-name parameter map; all required fixed and
/// estimated values must be present.
/// @param age Integer model age 0..15, with age 15 representing the plus group.
/// @return Mature female weight/fecundity contribution in the source
/// convention.
template <class T>
T mature_weight_at_age(const InputData &data, const ParameterValues<T> &p,
                       int age) {
  if (age < data.biology.first_mature_age)
    return T(0);
  const auto key =
      age_length_key(data, p, 0, age, maturity_integration_time(data));
  T sum = T(0);
  for (size_t b = 0; b < key.size(); ++b) {
    const double length = length_midpoint(data, b);
    const T weight = weight_at_length(p, 0, length);
    const T fecundity = weight * (p.at("Eggs/kg_inter_Fem_GP_1") +
                                  p.at("Eggs/kg_slope_wt_Fem_GP_1") * weight);
    sum += key[b] * maturity_at_length(p, length) * fecundity;
  }
  return sum;
}
/// Build a discrete block-activity cache key for a year; use only within one
/// evaluation to avoid reusing AD tape vertices.
/// @param data Validated original-input schema for the supported swordfish
/// assessment.
/// @param year Calendar year used to select observations, deviations or active
/// control blocks.
/// @return The typed result described above, with the module units and indexing
/// conventions.
inline std::vector<int> selectivity_period(const InputData &data, int year) {
  std::vector<int> active;
  for (const auto &ranges : data.blocks)
    for (size_t b = 0; b < ranges.size() / 2; ++b)
      active.push_back(year >= ranges[2 * b] && year <= ranges[2 * b + 1]);
  return active;
}
/// Resolve a fleet selectivity parameter with year-specific replacement blocks
/// and canonical parameter mapping.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active
/// evaluation tape.
/// @param data Validated original-input schema for the supported swordfish
/// assessment.
/// @param p Native canonical-name parameter map; all required fixed and
/// estimated values must be present.
/// @param fleet One-based source fleet identifier, not a zero-based vector
/// offset.
/// @param index One-based selectivity-parameter index in the canonical fleet
/// parameter name.
/// @param year Calendar year used to select observations, deviations or active
/// control blocks.
/// @return The typed result described above, with the module units and indexing
/// conventions.
template <class T>
T effective_selectivity_parameter(const InputData &data,
                                  const ParameterValues<T> &p, int fleet,
                                  int index, int year) {
  using std::exp;
  const std::string name =
      "SizeSel_fleet_" + std::to_string(fleet) + "_p" + std::to_string(index);
  const auto spec =
      std::find_if(data.parameters.begin(), data.parameters.end(),
                   [&](const Parameter &x) { return x.name == name; });
  if (spec == data.parameters.end())
    throw std::invalid_argument("Unknown selectivity parameter");
  T result = p.at(name);
  const int block = int(spec->columns.at(12));
  if (!block)
    return result;
  const auto &ranges = data.blocks.at(block - 1);
  for (size_t b = 0; b < ranges.size() / 2; ++b)
    if (year >= ranges[2 * b] && year <= ranges[2 * b + 1]) {
      const T adjustment = p.at(name + "_block_" + std::to_string(block) + "_" +
                                std::to_string(b + 1));
      switch (int(spec->columns.at(13))) {
      case 0:
        result *= exp(adjustment);
        break;
      case 1:
        result += adjustment;
        break;
      case 2:
        result = adjustment;
        break;
      default:
        throw std::invalid_argument("Unsupported selectivity block function");
      }
    }
  return result;
}
/// Build fleet length-selectivity curves for supported constant, logistic,
/// double-normal and mirror patterns, including active blocks.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active
/// evaluation tape.
/// @param data Validated original-input schema for the supported swordfish
/// assessment.
/// @param p Native canonical-name parameter map; all required fixed and
/// estimated values must be present.
/// @param year Calendar year used to select observations, deviations or active
/// control blocks.
/// @return The typed result described above, with the module units and indexing
/// conventions.
template <class T>
std::vector<std::vector<T>> length_selectivity(const InputData &data,
                                               const ParameterValues<T> &p,
                                               int year) {
  using std::exp;
  const size_t nbins = data.population_length_bins.size();
  std::vector<std::vector<T>> curves(data.fleets.size(),
                                     std::vector<T>(nbins, T(1)));
  for (size_t f = 0; f < curves.size(); ++f) {
    const auto &pattern = data.size_patterns[f];
    if (pattern[0] == 15) {
      if (pattern[3] < 1 || pattern[3] > int(f))
        throw std::invalid_argument("Selectivity mirror order");
      curves[f] = curves[pattern[3] - 1];
      continue;
    }
    if (pattern[0] == 0)
      continue;
    if (pattern[0] != 1 && pattern[0] != 24)
      throw std::invalid_argument("Selectivity pattern");
    const int n = pattern[0] == 1 ? 2 : 6;
    std::vector<T> sp;
    for (int i = 1; i <= n; ++i)
      sp.push_back(
          effective_selectivity_parameter(data, p, int(f) + 1, i, year));
    if (pattern[0] == 24 &&
        (sp[4] < -999 || sp[4] > -999 || sp[5] < -999 || sp[5] > -999))
      throw std::invalid_argument("This pattern-24 implementation requires "
                                  "supplied -999 endpoint settings");
    for (size_t b = 0; b < nbins; ++b) {
      const double length = length_midpoint(data, b);
      if (pattern[0] == 1) {
        if (!(sp[1] > 0))
          throw std::invalid_argument("Selectivity logistic width");
        curves[f][b] = logistic(T(std::log(19)) * (T(length) - sp[0]) / sp[1]);
      } else {
        const double width =
            data.population_length_bins[1] - data.population_length_bins[0];
        const T peak = sp[0],
                peak2 = peak + width +
                        (T(0.99 * length_midpoint(data, nbins - 1)) - peak -
                         width) *
                            logistic(sp[1]);
        const T t1 = T(length) - peak, t2 = T(length) - peak2;
        const T j1 = logistic(T(20) * t1 / (T(1) + (t1 < 0 ? -t1 : t1)));
        const T j2 = logistic(T(20) * t2 / (T(1) + (t2 < 0 ? -t2 : t2)));
        const T asc = exp(-t1 * t1 / exp(sp[2])),
                desc = exp(-t2 * t2 / exp(sp[3]));
        curves[f][b] = asc * (T(1) - j1) + j1 * (T(1) - j2 + desc * j2);
      }
    }
  }
  return curves;
}
/// Integrated selectivity in numbers and biomass, with selected mean length and
/// weight for a sex/age cell.
/// @details numbers is integrated number selectivity; biomass integrates
/// selectivity times individual kg weight. mean_length and mean_weight describe
/// selected fish, not the entire unselected population.
template <class T> struct SelectedAge {
  T numbers = T(0), biomass = T(0), mean_length = T(0), mean_weight = T(0);
};
/// Integrate fleet length selectivity and weight over an age-length key to
/// obtain selected-at-age summary quantities.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active
/// evaluation tape.
/// @param data Validated original-input schema for the supported swordfish
/// assessment.
/// @param p Native canonical-name parameter map; all required fixed and
/// estimated values must be present.
/// @param curve One fleet length-selectivity vector on the population
/// integration grid.
/// @param sex Zero-based sex index: 0 female, 1 male.
/// @param age Integer model age 0..15, with age 15 representing the plus group.
/// @param time Within-year time in years at which growth/length integration is
/// evaluated.
/// @return The typed result described above, with the module units and indexing
/// conventions.
template <class T>
SelectedAge<T>
selected_at_age(const InputData &data, const ParameterValues<T> &p,
                const std::vector<T> &curve, int sex, int age, double time) {
  const auto key = age_length_key(data, p, sex, age, time);
  if (curve.size() != key.size())
    throw std::invalid_argument("Selectivity grid dimension");
  SelectedAge<T> out;
  for (size_t b = 0; b < key.size(); ++b) {
    const T selected = key[b] * curve[b];
    const double l = length_midpoint(data, b);
    out.numbers += selected;
    out.biomass += selected * weight_at_length(p, sex, l);
    out.mean_length += selected * l;
  }
  if (!(out.numbers > 0))
    throw std::invalid_argument("Zero selected probability");
  out.mean_weight = out.biomass / out.numbers;
  out.mean_length /= out.numbers;
  return out;
}
} // namespace swordfish
