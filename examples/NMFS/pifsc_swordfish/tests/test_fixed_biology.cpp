#include "core/autodiff/model_gradient.hpp"
#include "swordfish/model.hpp"
#include "swordfish/reference_io.hpp"
#include "swordfish_input.hpp"
#include <iostream>
DECLARE_ADGRAPH();

/// Compare the entire objective/gradient, keeping the original equations as the
/// oracle.
void compare_paths(swordfish::SwordfishModel &model,
                   const std::vector<double> &point) {
  model.use_fixed_biology_cache = false;
  const auto original = quadra::evaluate_gradient(model, point);
  model.use_fixed_biology_cache = true;
  const auto cached = quadra::evaluate_gradient(model, point);
  if (std::abs(original.objective_value_m - cached.objective_value_m) > 1e-9)
    throw std::runtime_error("Fixed biology changed objective");
  for (size_t i = 0; i < point.size(); ++i)
    if (std::abs(original.gradient_m[i] - cached.gradient_m[i]) >
        1e-8 * (1 + std::abs(original.gradient_m[i])))
      throw std::runtime_error("Fixed biology changed gradient: " +
                               model.free[i].name);
}

/// Check cached/uncached agreement, cache invalidation and estimated-growth
/// fallback.
int main() {
  try {
    auto data = imported_swordfish_data();
    const auto native = swordfish::read_reference_parameters(
        data, "build/biology_fitted_parameters.txt",
        "build/population_fitted_deviations.txt");
    swordfish::SwordfishModel model(data, native.first, native.second);
    swordfish::FixedBiology fixed(model.data, model.baseline);
    if (!fixed.enabled())
      throw std::runtime_error(
          "Fixed biological parameters did not enable workspace");
    const auto point = model.initial_vector();
    compare_paths(model, point);
    auto moved = point;
    moved[0] += .01;
    compare_paths(model, moved);
    // A new fixed value must be recomputed rather than borrowed from a prior
    // call.
    model.baseline.at("L_at_Amin_Fem_GP_1") += .01;
    compare_paths(model, point);
    for (auto &spec : data.parameters)
      if (spec.name == "VonBert_K_Fem_GP_1")
        spec.columns[6] = 1;
    swordfish::SwordfishModel growth(data, native.first, native.second);
    swordfish::FixedBiology variable(growth.data, growth.baseline);
    if (variable.enabled())
      throw std::runtime_error("Estimated growth entered numeric-only cache");
    const auto start = growth.initial_vector();
    compare_paths(growth, start);
    const auto gradient = quadra::evaluate_gradient(growth, start);
    for (size_t i = 0; i < growth.free.size(); ++i)
      if (growth.free[i].name == "VonBert_K_Fem_GP_1") {
        auto plus = start, minus = start;
        plus[i] += 1e-6;
        minus[i] -= 1e-6;
        const double difference = (growth(plus) - growth(minus)) / 2e-6;
        if (std::abs(difference - gradient.gradient_m[i]) >
            1e-5 * (1 + std::abs(difference)))
          throw std::runtime_error("Estimated growth lost its derivative");
      }
    std::cout << "Fixed biology objective/all gradients, changed inputs and "
                 "estimated-growth fallback: PASS\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
