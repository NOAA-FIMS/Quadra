#include "core/autodiff/model_gradient.hpp"
#include "core/model/quadra_model.hpp"
#include "swordfish/population.hpp"
#include "swordfish_biology_input.hpp"
#include <iostream>

DECLARE_ADGRAPH();
using namespace swordfish;
/// Test helper for test population: require. See the typed signature and
/// assertions in the linked test source.
void require(bool ok, const char *message) {
  if (!ok)
    throw std::runtime_error(message);
}
/// Test helper for test population: near. See the typed signature and
/// assertions in the linked test source.
void near(double a, double b, double tol, const char *message) {
  require(std::isfinite(a) && std::isfinite(b) && std::abs(a - b) <= tol,
          message);
}
template <class T> SeasonalBiology<T> synthetic_biology() {
  /// Test helper for test population: synthetic biology. See the typed
  /// signature and assertions in the linked test source.
  SeasonalBiology<T> b;
  b.selectivity.resize(nfishing);
  b.selected_biomass.resize(nfishing);
  for (int f = 0; f < nfishing; ++f) {
    b.selectivity[f].fill(T(0));
    b.selected_biomass[f].fill(T(0));
  }
  b.selectivity[0][cell(0, 3)] = T(0.8);
  b.selected_biomass[0][cell(0, 3)] = T(8);
  b.selectivity[1][cell(0, 3)] = T(0.6);
  b.selected_biomass[1][cell(0, 3)] = T(12);
  return b;
}
/// Test helper for test population: PopulationDerivative. See the typed
/// signature and assertions in the linked test source.
class PopulationDerivative : public quadra::QuadraModel<PopulationDerivative> {
public:
  InputData data = imported_swordfish_biology_data();
  std::array<SeasonalBiology<double>, 4> fixed_biology;
  /// Test helper for test population: PopulationDerivative /
  /// PopulationDerivative. See the typed signature and assertions in the linked
  /// test source.
  PopulationDerivative() {
    const auto p = initial_values<double>(data);
    const auto curves = length_selectivity(data, p, data.dimensions.first_year);
    for (int s = 0; s < 4; ++s)
      fixed_biology[s] = seasonal_biology(data, p, curves, s + 1);
  }
  /// Test helper for test population: PopulationDerivative / parameter names
  /// impl. See the typed signature and assertions in the linked test source.
  std::vector<std::string> parameter_names_impl() const {
    return {"log_abundance", "initial_f", "early_deviation"};
  }
  /// Test helper for test population: PopulationDerivative / evaluate impl. See
  /// the typed signature and assertions in the linked test source.
  template <class T>
  T evaluate_impl(const std::vector<T> &x, quadra::ModelReportContext &) const {
    using std::exp;
    State<T> numbers = zero_state<T>();
    numbers[cell(0, 3)] = exp(x[0]);
    std::array<double, nfishing> obs{};
    obs[0] = 15;
    obs[1] = 120;
    const auto hybrid =
        hybrid_catch(data, numbers, synthetic_biology<T>(), obs, 0.25);
    auto p = initial_values<T>(data);
    p.at("InitF_seas_1_flt_6") = x[1];
    // These biological quantities do not depend on the tested parameters.
    // Their growth/selectivity derivatives are tested in test_biology.cpp.
    std::array<SeasonalBiology<T>, 4> biology;
    for (int s = 0; s < 4; ++s) {
      auto &dest = biology[s];
      const auto &source = fixed_biology[s];
      dest.selectivity.resize(nfishing);
      dest.selected_biomass.resize(nfishing);
      for (int j = 0; j < 32; ++j) {
        dest.weight_begin[j] = T(source.weight_begin[j]);
        for (int f = 0; f < nfishing; ++f) {
          dest.selectivity[f][j] = T(source.selectivity[f][j]);
          dest.selected_biomass[f][j] = T(source.selected_biomass[f][j]);
        }
      }
      for (int a = 0; a < nage; ++a)
        dest.mature_weight[a] = T(source.mature_weight[a]);
    }
    Deviations<T> dev;
    dev.emplace(1974, x[2]);
    const auto initial = initialize_population(data, p, dev, biology);
    return hybrid.apical_f[0] + T(0.5) * hybrid.apical_f[1] +
           T(0.001) * initial.january[cell(0, 1)] +
           T(0.00001) * initial.initial.ssb;
  }
};
/// Run test population; return zero on success and nonzero on a failed check or
/// invalid input.
int main() {
  try {
    auto data = imported_swordfish_biology_data();
    const auto b = synthetic_biology<double>();
    State<double> numbers{};
    numbers[cell(0, 3)] = 1000;
    std::array<double, nfishing> zero{};
    const auto none = hybrid_catch(data, numbers, b, zero, 0.25);
    for (int f = 0; f < nfishing; ++f) {
      near(none.apical_f[f], 0, 0, "Zero catch F");
      near(none.predicted_catch[f], 0, 0, "Zero catch prediction");
    }
    const double exposure_value = exposure(
        data.biology.natural_mortality[cell(0, 3)] + 0.8 * 0.15 + 0.6 * 0.25,
        0.25);
    std::array<double, nfishing> obs{};
    obs[0] = 1000 * 0.8 * 0.15 * exposure_value;
    obs[1] = 1000 * 12 * 0.25 * exposure_value;
    const auto fit = hybrid_catch(data, numbers, b, obs, 0.25);
    near(fit.apical_f[0], 0.15, 1e-5, "Number-catch F recovery");
    near(fit.apical_f[1], 0.25, 1e-5, "Biomass-catch F recovery");
    for (int f = 0; f < 2; ++f)
      near(fit.predicted_catch[f], obs[f], obs[f] * 1e-5,
           "Mixed-unit catches recovered");
    obs[0] = 1e9;
    obs[1] = 0;
    const auto excessive = hybrid_catch(data, numbers, b, obs, 0.25);
    require(std::isfinite(excessive.apical_f[0]) &&
                excessive.apical_f[0] <= data.max_f + 0.01,
            "Infeasible catches retain finite capped F");
    require(excessive.predicted_catch[0] < obs[0],
            "Infeasible catch must remain visibly unmatched");
    const auto empty = hybrid_catch(data, State<double>{}, b, obs, 0.25);
    near(empty.predicted_catch[0], 0, 0, "No fish cannot create catch");
    auto invalid = zero;
    invalid[0] = -1;
    bool rejected = false;
    try {
      hybrid_catch(data, numbers, b, invalid, 0.25);
    } catch (const std::invalid_argument &) {
      rejected = true;
    }
    require(rejected, "Negative catch rejection");
    const auto &a = data.recruitment_advanced;
    require(recruitment_bias(data, 1965) == 0, "Early ramp zero");
    near(recruitment_bias(data, 2000), a[8], 1e-14, "Full bias interval");
    require(recruitment_bias(data, 1974) > 0 &&
                recruitment_bias(data, 1974) < a[8],
            "Early partial bias");
    require(recruitment_bias(data, 2021) > 0 &&
                recruitment_bias(data, 2021) < a[8],
            "Recent partial bias");
    PopulationDerivative model;
    std::vector<double> x = {std::log(1000.), 0.7, 0.02};
    const auto ad = quadra::evaluate_gradient(model, x);
    quadra::ModelReportContext ctx;
    for (size_t j = 0; j < x.size(); ++j) {
      auto hi = x, lo = x;
      double step = 1e-5;
      hi[j] += step;
      lo[j] -= step;
      const double fd =
          (model.evaluate<double>(hi, ctx) - model.evaluate<double>(lo, ctx)) /
          (2 * step);
      near(ad.gradient_m[j], fd, 1e-6,
           "Hybrid/initialization automatic derivative");
    }
    std::cout
        << "Hybrid mixed-unit inversion, zero/infeasible catches, recruitment "
           "ramp and initialization/solver derivatives: PASS\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
