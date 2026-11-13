#include "core/autodiff/model_gradient.hpp"
#include "core/model/quadra_model.hpp"
#include "swordfish/likelihood.hpp"
#include "swordfish_biology_input.hpp"
#include <iostream>

DECLARE_ADGRAPH();
using namespace swordfish;
/// Test helper for test likelihood: require. See the typed signature and
/// assertions in the linked test source.
void require(bool ok, const char *m) {
  if (!ok)
    throw std::runtime_error(m);
}
/// Test helper for test likelihood: near. See the typed signature and
/// assertions in the linked test source.
void near(double a, double b, double tol, const char *m) {
  require(std::isfinite(a) && std::isfinite(b) && std::abs(a - b) <= tol, m);
}
/// Test helper for test likelihood: Terms. See the typed signature and
/// assertions in the linked test source.
class Terms : public quadra::QuadraModel<Terms> {
public:
  InputData data = imported_swordfish_biology_data();
  /// Test helper for test likelihood: Terms / parameter names impl. See the
  /// typed signature and assertions in the linked test source.
  std::vector<std::string> parameter_names_impl() const {
    return {"logit_probability", "log_sigma", "selectivity_peak",
            "common_log_scale"};
  }
  /// Test helper for test likelihood: Terms / evaluate impl. See the typed
  /// signature and assertions in the linked test source.
  template <class T>
  T evaluate_impl(const std::vector<T> &x, quadra::ModelReportContext &) const {
    using std::exp;
    using std::log;
    const T probability = logistic(x[0]);
    T nll = multinomial_deviance(
        std::vector<double>{0.25, 0.75},
        std::vector<T>{probability, T(1) - probability}, 20);
    auto p = initial_values<T>(data);
    p.at("SR_sigmaR") = exp(x[1]);
    p.at("SizeSel_fleet_1_p1") = x[2];
    Deviations<T> dev;
    dev.emplace(1975, T(0.2));
    dev.emplace(1990, T(-0.1));
    nll += recruitment_nll(data, p, dev) + selectivity_soft_bounds(data, p);
    const std::vector<T> log_base = {T(std::log(2.)) + x[3],
                                     T(std::log(3.)) + x[3]};
    const std::vector<double> obs = {1, 2}, se = {0.2, 0.4};
    const T q = profiled_log_q(log_base, obs, se, {true, true}, false);
    nll += lognormal_index_nll(log_base[0] + q, obs[0], se[0], 1);
    nll += lognormal_index_nll(log_base[1] + q, obs[1], se[1], 1);
    return nll;
  }
};
/// Run test likelihood; return zero on success and nonzero on a failed check or
/// invalid input.
int main() {
  try {
    auto data = imported_swordfish_biology_data();
    near(likelihood_weight(data, 1, 24, 5), 0, 0,
         "Survey24 excluded by lambda");
    near(likelihood_weight(data, 1, 27, 5), 0, 0,
         "Survey27 excluded by lambda");
    near(likelihood_weight(data, 4, 1, 5), 0.5, 0,
         "Length lambda propagated to later phases");
    near(likelihood_weight(data, 9, 6, 5), 0, 0, "Equilibrium catch excluded");
    near(variance_adjustment(data, 1, 24, 0), 0.11, 0, "Index SE adjustment");
    const std::vector<double> probability = {0.25, 0.75};
    near(multinomial_deviance(probability, probability, 20.), 0, 1e-14,
         "Perfect-fit multinomial offset");
    require(multinomial_deviance(probability, std::vector<double>{0.5, 0.5},
                                 20.) > 0,
            "Imperfect multinomial fit");
    near(stabilized_catch_nll(100., 100., 0.05), 0, 1e-20,
         "Perfect-fit catch objective");
    const auto q =
        profiled_log_q(std::vector<double>{std::log(2.), std::log(3.)}, {1, 2},
                       {0.2, 0.4}, {true, true}, false);
    const double expected = (std::log(0.5) / 0.04 + std::log(2. / 3) / 0.16) /
                            (1 / 0.04 + 1 / 0.16);
    near(q, expected, 1e-14, "Inverse-variance Q profile");
    near(profiled_log_q(std::vector<double>{0., 0.}, {1, 1}, {0.2, 0.2},
                        {false, false}, false),
         0, 0, "No used Q observations");
    Terms model;
    std::vector<double> x = {-0.4, std::log(0.42), 147.296, 0.3};
    const auto gradient = quadra::evaluate_gradient(model, x);
    quadra::ModelReportContext context;
    for (size_t j = 0; j < x.size(); ++j) {
      auto hi = x, lo = x;
      const double step = 1e-5;
      hi[j] += step;
      lo[j] -= step;
      const double fd = (model.evaluate<double>(hi, context) -
                         model.evaluate<double>(lo, context)) /
                        (2 * step);
      near(gradient.gradient_m[j], fd, 1e-7, "Likelihood term AD derivative");
    }
    near(gradient.gradient_m[0], 20 * (logistic(x[0]) - 0.25), 1e-12,
         "Analytic multinomial derivative");
    double bias = 0;
    for (int y = 1965; y <= 2021; ++y)
      bias += recruitment_bias(data, y);
    near(gradient.gradient_m[1], bias - 0.05 / (0.42 * 0.42), 1e-12,
         "Analytic recruitment sigma derivative");
    near(gradient.gradient_m[3], 0, 1e-12,
         "Profiled Q removes common scale from index objective");
    std::cout << "Likelihood offsets, phase weights, Q profiling and "
                 "analytic/Quadra derivatives: PASS\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
