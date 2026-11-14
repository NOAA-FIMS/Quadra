/** @file fixed_biology.hpp
 * Ordinary numeric biology shared by population and observations in one
 * evaluation. Never retains AD scalars or crosses a tape boundary. Estimated
 * biology uses the original differentiable calculations instead.
 */
#pragma once
#include "biology.hpp"

namespace swordfish {
/// Numeric sex/age/length probabilities for one within-year sampling time.
using NumericLengthKeys = std::array<std::vector<std::vector<double>>, 2>;

/// A small evaluation-scoped workspace for biology whose parameters are fixed.
/// Keeping the dependency check explicit makes the fast path conservative:
/// estimating any biological parameter disables the entire numeric workspace.
/// The small native parameter snapshot is owned; input data must outlive this
/// workspace. Values are recomputed for each evaluation, so changed inputs
/// cannot leave stale persistent state. Only ordinary doubles are stored.
class FixedBiology {
  const InputData &data;
  ParameterValues<double> parameters;
  bool usable;
  std::map<double, NumericLengthKeys> keys;
  std::map<double, State<double>> beginning_weights;

public:
  std::array<std::vector<double>, 2> length_weights;
  std::array<double, nage> mature_weights{};
  /// Check biological estimation flags before calculating any cached quantity.
  FixedBiology(const InputData &d, const ParameterValues<double> &p,
               bool requested = true)
      : data(d), parameters(p), usable(requested) {
    for (const auto &spec : data.parameters) {
      const bool biological = spec.name.find("_GP_1") != std::string::npos ||
                              spec.name == "CohortGrowDev";
      if (biological && spec.columns[6] > 0)
        usable = false;
    }
    if (!usable)
      return;
    for (int sex = 0; sex < 2; ++sex)
      for (size_t bin = 0; bin < data.population_length_bins.size(); ++bin)
        length_weights[sex].push_back(
            weight_at_length(parameters, sex, length_midpoint(data, bin)));
    for (int age = 0; age < nage; ++age)
      mature_weights[age] = mature_weight_at_age(data, parameters, age);
  }
  /// True only when all biological parameter dependencies are fixed.
  bool enabled() const { return usable; }
  /// Compute a time's numeric age-length keys once, shared across consumers.
  const NumericLengthKeys &length_keys(double time) {
    if (!usable)
      throw std::logic_error("Fixed biology requested for estimated biology");
    auto found = keys.find(time);
    if (found == keys.end()) {
      NumericLengthKeys value;
      for (int sex = 0; sex < 2; ++sex)
        for (int age = 0; age < nage; ++age)
          value[sex].push_back(
              age_length_key(data, parameters, sex, age, time));
      found = keys.emplace(time, std::move(value)).first;
    }
    return found->second;
  }
  /// Compute beginning-season individual weights once for each requested time.
  const State<double> &weights_at_age(double time) {
    if (!usable)
      throw std::logic_error("Fixed weights requested for estimated biology");
    auto found = beginning_weights.find(time);
    if (found == beginning_weights.end()) {
      auto value = zero_state<double>();
      for (int sex = 0; sex < 2; ++sex)
        for (int age = 0; age < nage; ++age)
          value[cell(sex, age)] =
              weight_at_age(data, parameters, sex, age, time);
      found = beginning_weights.emplace(time, std::move(value)).first;
    }
    return found->second;
  }
};
} // namespace swordfish
