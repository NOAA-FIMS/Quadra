/**
 * @file likelihood.hpp
 * @brief SS3-compatible prediction and joint likelihood accounting.
 *
 * Observation units follow their fleet definitions. Log catchability is
 * profiled analytically; native parameter maps may therefore contain Q values
 * that are not independent fitting coordinates. Component functions return
 * unweighted contributions unless explicitly described otherwise;
 * evaluate_likelihood applies the phase/fleet lambdas and keeps excluded length
 * samples out of the target.
 */
#pragma once
#include "population.hpp"

namespace swordfish {
/// Resolve a component/fleet likelihood lambda at the requested estimation
/// phase.
/// @param data Validated original-input schema for the supported swordfish
/// assessment.
/// @param component SS3 likelihood-component code, as stored in lambda_changes.
/// @param fleet One-based source fleet identifier, not a zero-based vector
/// offset.
/// @param phase Estimation phase used to resolve persistent lambda changes.
/// @return Effective lambda after applicable changes persist through the
/// selected phase.
inline double likelihood_weight(const InputData &data, int component, int fleet,
                                int phase) {
  const int used_phase = std::min(phase, data.max_lambda_phase);
  double weight = 1;
  // Input changes persist through later phases, in the supplied row order.
  for (const auto &row : data.lambda_changes)
    if (int(row[0]) == component && int(row[1]) == fleet &&
        int(row[2]) <= used_phase)
      weight = row[3];
  return weight;
}
/// Apply a parsed variance adjustment for an observation kind and fleet to its
/// initial error value.
/// @param data Validated original-input schema for the supported swordfish
/// assessment.
/// @param kind SS3 variance-adjustment kind code.
/// @param fleet One-based source fleet identifier, not a zero-based vector
/// offset.
/// @param initial Default adjustment/error value when no matching control row
/// exists.
/// @return Configured replacement value or initial when no row matches.
inline double variance_adjustment(const InputData &data, int kind, int fleet,
                                  double initial) {
  double value = initial;
  for (const auto &row : data.variance_adjustments)
    if (int(row[0]) == kind && int(row[1]) == fleet)
      value = row[2];
  return value;
}
/// Return a checked reference to one historical year-season row; reject missing
/// or out-of-range calendar positions.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active
/// evaluation tape.
/// @param run Completed historical population trajectory; references into it
/// must not outlive it.
/// @param data Validated original-input schema for the supported swordfish
/// assessment.
/// @param year Calendar year used to select observations, deviations or active
/// control blocks.
/// @param season One-based quarter, 1..4.
/// @return Borrowed reference into run; destruction/mutation of run invalidates
/// the reference.
template <class T>
const PopulationSeason<T> &population_row(const PopulationRun<T> &run,
                                          const InputData &data, int year,
                                          int season) {
  if (year < data.dimensions.first_year || year > data.dimensions.last_year ||
      season < 1 || season > 4)
    throw std::invalid_argument("Observation outside population run");
  return run.seasons.at(
      size_t((year - data.dimensions.first_year) * 4 + season - 1));
}
/// Resolved season, within-year observation time, age-length-key time and
/// elapsed seasonal survival duration.
/// @details season is one-based. time and alk_time are calendar-year fractions;
/// elapsed is time since the season start. Negative elapsed is an intentional
/// sentinel for season-integrated sampling, not a negative survival duration.
struct ObservationTime {
  int season;
  double time, alk_time, elapsed;
};
/// Translate a source observation month and fleet timing setting into the
/// season and exposure times used for index and composition prediction.
/// @param data Validated original-input schema for the supported swordfish
/// assessment.
/// @param fleet One-based source fleet identifier, not a zero-based vector
/// offset.
/// @param month Source observation month, in [1,13); fractional months retain
/// their supplied timing.
/// @return The typed result described above, with the module units and indexing
/// conventions.
inline ObservationTime observation_time(const InputData &data, int fleet,
                                        double month) {
  if (fleet < 1 || fleet > 27 || !(month >= 1 && month < 13))
    throw std::invalid_argument("Observation timing");
  const double time = (month - 1) / 12;
  const int season = int(std::floor(time * 4 + 1e-10)) + 1;
  const double start = (season - 1) * 0.25;
  if (data.fleets.at(fleet - 1).timing == -1)
    return {season, start + 0.125, start + 0.125, -1};
  return {season, time, std::floor(time / 0.125 + 1e-10) * 0.125, time - start};
}
/// Sex-specific cached age-length keys for a single within-year time in one
/// likelihood evaluation.
/// @details Keys are indexed first by sex, then age, then length bin. A cached
/// key must have the same parameters and observation integration time as its
/// consumer; AD keys cannot cross an evaluation/tape boundary.
template <class T>
using LengthKeys = std::array<std::vector<std::vector<T>>, 2>;
/// Project selected surviving numbers into sex-specific sampled length
/// frequencies for an observation time.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active
/// evaluation tape.
/// @param data Validated original-input schema for the supported swordfish
/// assessment.
/// @param p Native canonical-name parameter map; all required fixed and
/// estimated values must be present.
/// @param row One historical population-season row with beginning numbers and
/// fitted fishing rates.
/// @param curve One fleet length-selectivity vector on the population
/// integration grid.
/// @param time Resolved observation timing; negative elapsed selects the
/// supported season-integrated convention.
/// @param cached Optional age-length keys from this same parameter evaluation
/// and time; nullptr recomputes them.
/// @param fixed Optional ordinary numeric keys when all biological dependencies
/// are fixed; these never retain AD values across evaluations.
/// @return The typed result described above, with the module units and indexing
/// conventions.
template <class T>
std::array<std::vector<T>, 2>
sampled_lengths(const InputData &data, const ParameterValues<T> &p,
                const PopulationSeason<T> &row, const std::vector<T> &curve,
                const ObservationTime &time,
                const LengthKeys<T> *cached = nullptr,
                const NumericLengthKeys *fixed = nullptr) {
  using std::exp;
  std::array<std::vector<T>, 2> length;
  for (auto &v : length)
    v.assign(data.population_length_bins.size(), T(0));
  for (int sex = 0; sex < 2; ++sex)
    for (int a = 0; a < nage; ++a) {
      const int j = cell(sex, a);
      T z = T(data.biology.natural_mortality[j]);
      for (const auto &f : row.hybrid.fishing)
        z += f[j];
      const T available =
          row.beginning[j] *
          (time.elapsed < 0 ? exposure(z, 0.25) : exp(-z * time.elapsed));
      auto accumulate = [&](const auto &key) {
        using KeyScalar = typename std::decay_t<decltype(key)>::value_type;
        for (size_t b = 0; b < key.size(); ++b) {
          if constexpr (std::is_same_v<KeyScalar, double>)
            if (key[b] == 0.)
              continue;
          length[sex][b] += available * key[b] * curve[b];
        }
      };
      if (fixed)
        accumulate(fixed->at(sex).at(a));
      else if (cached)
        accumulate(cached->at(sex).at(a));
      else
        accumulate(age_length_key(data, p, sex, a, time.alk_time));
    }
  return length;
}
/// Analytically profile log catchability as an inverse-variance-weighted log
/// residual mean, using the exact source stabilization.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active
/// evaluation tape.
/// @param log_base Log model predictions before catchability, already using the
/// required stabilizer.
/// @param observed Strictly positive index observations aligned with log_base.
/// @param se Positive log-scale standard errors; inverse squared values define
/// the profiling weights.
/// @param use Inclusion flags aligned with observations and error vectors.
/// @param bias_adjust Whether to apply the supported lognormal catchability
/// bias correction.
/// @return Analytical optimum log catchability; zero when no included
/// observations exist.
template <class T>
T profiled_log_q(const std::vector<T> &log_base,
                 const std::vector<double> &observed,
                 const std::vector<double> &se, const std::vector<bool> &use,
                 bool bias_adjust) {
  using std::log;
  if (log_base.size() != observed.size() || se.size() != observed.size() ||
      use.size() != se.size())
    throw std::invalid_argument("Profiled Q dimensions");
  // Profiling Q at every evaluation retains its dependence on the population
  // parameters in AD. Treating saved report Q as fixed would change gradients.
  T sum = T(0);
  double precision = 0;
  size_t n = 0;
  for (size_t i = 0; i < observed.size(); ++i)
    if (use[i]) {
      if (!(observed[i] > 0 && se[i] > 0))
        throw std::invalid_argument("Profiled Q observation");
      const double w = 1 / (se[i] * se[i]);
      sum += (T(std::log(observed[i])) - log_base[i]) * w;
      precision += w;
      ++n;
    }
  if (n == 0)
    return T(0);
  if (bias_adjust)
    sum += T(0.5 * n);
  return sum / precision;
}
/// Evaluate the index lognormal negative log-likelihood under the configured
/// standard-deviation offset convention.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active
/// evaluation tape.
/// @param log_expected Expected observation on the natural-log scale.
/// @param observed Strictly positive index observation on its native scale.
/// @param se Positive log-scale standard error for this observation.
/// @param sd_offset Configured multiplier of the log-standard-error term in the
/// likelihood.
/// @return Unweighted lognormal contribution, including the configured SD
/// offset.
template <class T>
T lognormal_index_nll(const T &log_expected, double observed, double se,
                      int sd_offset) {
  using std::log;
  const T residual = (T(std::log(observed)) - log_expected) / se;
  return T(0.5) * residual * residual + T(sd_offset * std::log(se));
}
/// Calculate sample-size-weighted multinomial composition deviance against a
/// perfect-fit offset after probability stabilization.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active
/// evaluation tape.
/// @param observed Observed length-bin proportions after the supported
/// stabilization/normalization.
/// @param expected Strictly positive expected bin proportions aligned with
/// observed.
/// @param n Positive effective sample size after the supplied composition
/// weighting.
/// @return Unweighted multinomial perfect-fit-offset deviance.
template <class T>
T multinomial_deviance(const std::vector<double> &observed,
                       const std::vector<T> &expected, double n) {
  using std::log;
  if (observed.size() != expected.size() || !(n > 0))
    throw std::invalid_argument("Multinomial dimensions/sample size");
  T value = T(0);
  for (size_t b = 0; b < observed.size(); ++b)
    if (observed[b] > 0) {
      if (!(expected[b] > 0))
        throw std::invalid_argument("Zero expected composition probability");
      value +=
          T(n * observed[b]) * (T(std::log(observed[b])) - log(expected[b]));
    }
  return value;
}
/// Evaluate SS3's stabilized log catch residual penalty for observed and
/// predicted fleet-unit catch.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active
/// evaluation tape.
/// @param observed Strictly positive fleet-unit observed catch.
/// @param predicted Predicted catch in the same fleet units as observed.
/// @param se Positive log-scale catch standard error.
/// @return Unweighted stabilized log catch residual penalty.
template <class T>
T stabilized_catch_nll(double observed, const T &predicted, double se) {
  using std::log;
  if (!(observed > 0 && se > 0))
    throw std::invalid_argument("Catch likelihood inputs");
  const T residual =
      (T(std::log(1.1 * observed)) - log(predicted + T(0.1 * observed))) / se;
  return T(0.5) * residual * residual;
}
/// Evaluate the supported symmetric-beta soft-bound penalty on estimated
/// selectivity parameters.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active
/// evaluation tape.
/// @param data Validated original-input schema for the supported swordfish
/// assessment.
/// @param p Native canonical-name parameter map; all required fixed and
/// estimated values must be present.
/// @return The typed result described above, with the module units and indexing
/// conventions.
template <class T>
T selectivity_soft_bounds(const InputData &data, const ParameterValues<T> &p) {
  using std::log;
  T penalty = T(0);
  for (const auto &spec : data.parameters)
    if (spec.name.rfind("SizeSel_", 0) == 0 && spec.columns[6] > 0) {
      const double lo = spec.columns[0], hi = spec.columns[1];
      const T value = p.at(spec.name);
      const double mu = -0.001 * (std::log((hi - lo) * 0.5) + std::log(0.5));
      penalty -= T(mu) + T(0.001) * log(value - T(lo) + T(0.0001)) +
                 T(0.001) * log(T(1) - (value - T(lo) - T(0.0001)) / (hi - lo));
    }
  return penalty;
}
/// Sum squared historical log recruitment deviations and the bias-weighted
/// log-sigma contribution for the joint likelihood target.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active
/// evaluation tape.
/// @param data Validated original-input schema for the supported swordfish
/// assessment.
/// @param p Native canonical-name parameter map; all required fixed and
/// estimated values must be present.
/// @param dev Year-to-log-recruitment-deviation map; the caller enforces the
/// main-period constraint.
/// @return The typed result described above, with the module units and indexing
/// conventions.
template <class T>
T recruitment_nll(const InputData &data, const ParameterValues<T> &p,
                  const Deviations<T> &dev) {
  using std::log;
  // The supplied model has independent Gaussian deviations, no autocorrelation.
  if (p.at("SR_autocorr") < 0 || p.at("SR_autocorr") > 0)
    throw std::invalid_argument("Unsupported recruitment autocorrelation");
  const T sigma = p.at("SR_sigmaR");
  if (!(sigma > 0))
    throw std::invalid_argument("Recruitment sigma");
  T sum = T(0);
  double bias_count = 0;
  for (int year = int(data.recruitment_advanced[0]); year <= data.recdev_end;
       ++year) {
    const auto it = dev.find(year);
    const T d = it == dev.end() ? T(0) : it->second;
    sum += d * d / (T(2) * sigma * sigma);
    bias_count += recruitment_bias(data, year);
  }
  return sum + T(data.sd_offset * bias_count) * log(sigma);
}
/// Index observation prediction, profiled catchability, residual, error and
/// unweighted likelihood contribution.
/// @details observation/prediction retain index fleet units, se is log-scale,
/// and log_q is analytically profiled catchability. nll is an unweighted
/// contribution; component lambdas are applied in the aggregate objective.
template <class T> struct IndexFit {
  int year = 0, season = 0, fleet = 0;
  double observation = 0, se = 0;
  T base = T(0), log_q = T(0), prediction = T(0), deviation = T(0), nll = T(0);
};
/// Length-sample prediction and diagnostics, including effective sample size,
/// observed/predicted means and deviance.
/// @details observed/expected are stabilized normalized bin proportions.
/// sample_size is the supplied effective N, means are cm, and nll is the
/// unweighted multinomial perfect-fit-offset contribution.
template <class T> struct LengthFit {
  int year = 0, season = 0, fleet = 0;
  double sample_size = 0, observed_mean = 0;
  std::vector<double> observed;
  std::vector<T> expected;
  T expected_mean = T(0), nll = T(0);
};
/// Joint objective total, named weighted components, fitted indices and
/// included length-sample diagnostics.
/// @details total is the weighted joint objective. components keeps named terms
/// for SS3 comparisons, while indices/lengths retain per-observation
/// diagnostics. Profiled Q values do not add independent optimization
/// coordinates.
template <class T> struct LikelihoodResult {
  std::map<std::string, T> components;
  std::vector<IndexFit<T>> indices;
  std::vector<LengthFit<T>> lengths;
  T total = T(0);
};
/// Calculate catch, index, length-composition, recruitment and soft-bound
/// components from a population trajectory; preserve all configured lambdas.
/// @tparam T Ordinary or AD scalar; all returned AD values belong to the active
/// evaluation tape.
/// @param data Validated original-input schema for the supported swordfish
/// assessment.
/// @param p Native canonical-name parameter map; all required fixed and
/// estimated values must be present.
/// @param deviations Historical fitted log recruitment deviations, including
/// early age-structure effects.
/// @param run Completed historical population trajectory; references into it
/// must not outlive it.
/// @param phase Estimation phase used to resolve persistent lambda changes.
/// @param fixed Optional numeric biology from this same evaluation and inputs.
/// @return The typed result described above, with the module units and indexing
/// conventions.
template <class T>
LikelihoodResult<T>
evaluate_likelihood(const InputData &data, const ParameterValues<T> &p,
                    const Deviations<T> &deviations,
                    const PopulationRun<T> &run, int phase = 5,
                    FixedBiology *fixed = nullptr) {
  using std::exp;
  using std::log;
  LikelihoodResult<T> out;
  for (const auto &name :
       {"Catch", "Equil_catch", "Survey", "Length_comp", "Recruitment",
        "InitEQ_Regime", "Forecast_Recruitment", "Parm_priors",
        "Parm_softbounds", "Parm_devs", "Crash_Pen"})
    out.components.emplace(name, T(0));
  for (const auto &spec : data.parameters)
    if (spec.columns[5] != 0)
      throw std::invalid_argument("Unsupported nonzero parameter prior");
  std::map<std::vector<int>, std::vector<std::vector<T>>> curves;
  std::map<int, LengthKeys<T>> keys;
  auto curve_for = [&](int year, int fleet) -> const std::vector<T> & {
    const auto key = selectivity_period(data, year);
    if (!curves.count(key))
      curves.emplace(key, length_selectivity(data, p, year));
    return curves.at(key).at(fleet - 1);
  };
  auto sample_for = [&](const PopulationSeason<T> &row,
                        const std::vector<T> &curve,
                        const ObservationTime &time) {
    if (fixed && fixed->enabled())
      return sampled_lengths(data, p, row, curve, time,
                             static_cast<const LengthKeys<T> *>(nullptr),
                             &fixed->length_keys(time.alk_time));
    const int slot = int(std::round(time.alk_time / 0.125));
    if (!keys.count(slot)) {
      LengthKeys<T> value;
      for (int sex = 0; sex < 2; ++sex)
        for (int a = 0; a < nage; ++a)
          value[sex].push_back(age_length_key(data, p, sex, a, time.alk_time));
      keys.emplace(slot, std::move(value));
    }
    return sampled_lengths(data, p, row, curve, time, &keys.at(slot));
  };
  // Profile catchability separately for each survey, retaining observations in
  // input order.
  for (const auto &x : data.indices) {
    const int year = std::abs(x.year), f = x.fleet;
    const auto time = observation_time(data, f, x.month);
    const auto &row = population_row(run, data, year, time.season);
    const auto &setup = data.index_specs.at(f - 1);
    if (setup[2] != 0)
      throw std::invalid_argument(
          "Only supplied lognormal index errors supported");
    IndexFit<T> fit;
    fit.year = year;
    fit.season = time.season;
    fit.fleet = f;
    fit.observation = x.value;
    fit.se = x.log_se + variance_adjustment(data, 1, f, 0);
    if (fit.se <= 0)
      fit.se = 0.001;
    if (setup[1] == 33) {
      if (time.season == 1)
        fit.base =
            year == data.dimensions.first_year
                ? adjusted_recruitment(data, p, deviations, year - 1,
                                       run.initialization.initial.recruitment)
                : population_row(run, data, year - 1, 3).recruits;
      else
        fit.base = population_row(run, data, year, 3).recruits;
    } else {
      if (setup[1] != 0 && setup[1] != 1)
        throw std::invalid_argument("Unsupported index units");
      const auto length = sample_for(row, curve_for(year, f), time);
      for (int sex = 0; sex < 2; ++sex)
        for (size_t b = 0; b < length[sex].size(); ++b)
          fit.base += length[sex][b] *
                      (setup[1] == 1
                           ? weight_at_length(p, sex, length_midpoint(data, b))
                           : T(1));
    }
    out.indices.push_back(std::move(fit));
  }
  for (int f = 1; f <= 27; ++f) {
    std::vector<T> log_base;
    std::vector<double> observed, se;
    std::vector<bool> use;
    for (size_t i = 0; i < data.indices.size(); ++i)
      if (data.indices[i].fleet == f) {
        const auto &fit = out.indices[i];
        log_base.push_back(log(fit.base + T(0.000001)));
        observed.push_back(fit.observation);
        se.push_back(fit.se);
        use.push_back(data.indices[i].use);
      }
    if (log_base.empty())
      continue;
    const auto q =
        std::find_if(data.q_setup.begin(), data.q_setup.end(),
                     [&](const auto &row) { return int(row[0]) == f; });
    if (q == data.q_setup.end() || (*q)[1] != 1 || (*q)[3] != 0 || (*q)[5] != 1)
      throw std::invalid_argument("Unsupported catchability setup");
    const T log_q = profiled_log_q(log_base, observed, se, use, (*q)[4] != 0);
    for (size_t i = 0; i < data.indices.size(); ++i)
      if (data.indices[i].fleet == f) {
        auto &fit = out.indices[i];
        fit.log_q = log_q;
        const T log_predicted = log(fit.base + T(0.000001)) + log_q;
        fit.prediction = exp(log_predicted);
        fit.deviation = T(std::log(fit.observation)) - log_predicted;
        fit.nll = lognormal_index_nll(log_predicted, fit.observation, fit.se,
                                      data.sd_offset);
        if (data.indices[i].use)
          out.components.at("Survey") +=
              fit.nll * likelihood_weight(data, 1, f, phase);
      }
  }
  for (const auto &x : data.compositions) {
    if (!x.use)
      continue;
    const int f = x.fleet, year = std::abs(x.year);
    const auto &ctl = data.length_controls.at(f - 1);
    if (x.sex != 0 || x.partition != 0 || ctl[0] >= 0 || ctl[2] != 0 ||
        ctl[3] != 0 || ctl[4] != 0 || ctl[5] != 0)
      throw std::invalid_argument("Composition implementation requires "
                                  "supplied combined-sex multinomial setup");
    const auto time = observation_time(data, f, x.month);
    const auto sample = sample_for(population_row(run, data, year, time.season),
                                   curve_for(year, f), time);
    LengthFit<T> fit;
    fit.year = year;
    fit.season = time.season;
    fit.fleet = f;
    fit.sample_size =
        std::max(ctl[6], x.sample_size * variance_adjustment(data, 4, f, 1));
    double obs_sum = 0;
    T exp_sum = T(0);
    for (size_t b = 0; b < 51; ++b) {
      obs_sum += x.values[b];
      exp_sum += sample[0][b] + sample[1][b];
    }
    if (!(obs_sum > 0 && exp_sum > 0))
      throw std::invalid_argument("Empty composition");
    for (size_t b = 0; b < 51; ++b) {
      const double observed =
          (x.values[b] / obs_sum + ctl[1]) / (1 + 51 * ctl[1]);
      const T expected = ((sample[0][b] + sample[1][b]) / exp_sum + T(ctl[1])) /
                         (1 + 51 * ctl[1]);
      fit.observed.push_back(observed);
      fit.expected.push_back(expected);
      fit.observed_mean += observed * length_midpoint(data, b);
      fit.expected_mean += expected * length_midpoint(data, b);
    }
    fit.nll = multinomial_deviance(fit.observed, fit.expected, fit.sample_size);
    out.components.at("Length_comp") +=
        fit.nll * likelihood_weight(data, 4, f, phase);
    out.lengths.push_back(std::move(fit));
  }
  for (const auto &x : data.catches)
    if (x.year != -999 && x.value > 0) {
      const auto &row = population_row(run, data, x.year, x.season);
      out.components.at("Catch") +=
          stabilized_catch_nll(x.value, row.hybrid.predicted_catch[x.fleet - 1],
                               x.log_se) *
          likelihood_weight(data, 8, x.fleet, phase);
    }
  // Initial equilibrium catch has zero lambda for the sole positive
  // observation.
  for (const auto &x : data.catches)
    if (x.year == -999 && x.value > 0 &&
        likelihood_weight(data, 9, x.fleet, phase) != 0)
      throw std::invalid_argument(
          "Nonzero equilibrium catch likelihood not yet supported");
  out.components.at("Recruitment") = recruitment_nll(data, p, deviations) *
                                     likelihood_weight(data, 10, 1, phase);
  out.components.at("Parm_softbounds") = selectivity_soft_bounds(data, p);
  for (const auto &component : out.components)
    out.total += component.second;
  return out;
}
} // namespace swordfish
