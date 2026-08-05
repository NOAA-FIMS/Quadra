#include "../include/quadra/stats.hpp"
#include "../core/autodiff/model_gradient.hpp"
#include "../core/model/quadra_model.hpp"

#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

DECLARE_ADGRAPH();

namespace {

constexpr double tolerance = 1e-11;

bool close(double lhs, double rhs, double tol = tolerance) {
  return std::abs(lhs - rhs) <= tol * (1.0 + std::abs(rhs));
}

double manual_ar1_logpdf(const std::vector<double> &x, double mean,
                         double phi, double innovation_sd) {
  if (x.empty()) {
    return 0.0;
  }
  const double marginal_sd = innovation_sd / std::sqrt(1.0 - phi * phi);
  double ans = quadra::dnorm(x[0], mean, marginal_sd, true);
  for (std::size_t i = 1; i < x.size(); ++i) {
    ans += quadra::dnorm(x[i], mean + phi * (x[i - 1] - mean),
                        innovation_sd, true);
  }
  return ans;
}

struct Ar1PhiModel : public quadra::QuadraModel<Ar1PhiModel> {
  std::vector<double> observations;
  double mean;
  double innovation_sd;

  Ar1PhiModel(std::vector<double> observations_in, double mean_in,
              double innovation_sd_in)
      : observations(std::move(observations_in)), mean(mean_in),
        innovation_sd(innovation_sd_in) {}

  std::vector<std::string> parameter_names_impl() const {
    return {"unconstrained_phi"};
  }

  template <typename T>
  T evaluate_impl(const std::vector<T> &p,
                  quadra::ModelReportContext &) const {
    std::vector<T> x;
    x.reserve(observations.size());
    for (double value : observations) {
      x.push_back(T(value));
    }
    const T phi = quadra::stats::correlation_from_unconstrained(p[0]);
    return quadra::stats::ar1_stationary_nll(
        x, T(mean), phi, T(innovation_sd));
  }

  double nll(double unconstrained_phi) const {
    quadra::ModelReportContext context;
    return evaluate_impl(std::vector<double>{unconstrained_phi}, context);
  }
};

double ar1_nll_at(const Ar1PhiModel &model, double unconstrained_phi) {
  return model.nll(unconstrained_phi);
}

} // namespace

int main() {
  using quadra::stats::ar1_marginal_logpdf;
  using quadra::stats::ar1_stationary_logpdf;

  const double mean = 0.3;
  const double phi = 0.65;
  const double innovation_sd = 0.8;
  const double marginal_sd = innovation_sd / std::sqrt(1.0 - phi * phi);
  const std::vector<double> x = {-0.2, 0.4, 1.1, 0.7};

  const double expected = manual_ar1_logpdf(x, mean, phi, innovation_sd);
  assert(close(ar1_stationary_logpdf(x, mean, phi, innovation_sd), expected));
  assert(close(ar1_marginal_logpdf(x, mean, phi, marginal_sd), expected));

  const std::vector<double> empty;
  assert(ar1_stationary_logpdf(empty, mean, phi, innovation_sd) == 0.0);

  const double iid_expected =
      quadra::stats::iid_normal_logpdf(x, mean, innovation_sd);
  assert(close(ar1_stationary_logpdf(x, mean, 0.0, innovation_sd),
               iid_expected));

  const std::vector<double> innovations = {-0.3, 0.2, 0.8, -0.1, 0.4};
  const double theta = -0.45;
  const auto ma1 =
      quadra::stats::ma1_from_innovations(innovations, mean, theta);
  assert(ma1.size() == innovations.size() - 1);
  for (std::size_t i = 0; i < ma1.size(); ++i) {
    assert(close(ma1[i], mean + innovations[i + 1] + theta * innovations[i]));
  }
  assert(close(quadra::stats::ma1_innovations_logpdf(innovations,
                                                      innovation_sd),
               quadra::stats::iid_normal_logpdf(innovations, 0.0,
                                                innovation_sd)));

  Ar1PhiModel model{x, mean, innovation_sd};
  const double unconstrained_phi = 0.4;
  const auto gradient =
      quadra::evaluate_gradient(model, std::vector<double>{unconstrained_phi});
  const double step = 1e-6;
  const double finite_difference =
      (ar1_nll_at(model, unconstrained_phi + step) -
       ar1_nll_at(model, unconstrained_phi - step)) /
      (2.0 * step);
  assert(gradient.gradient_m.size() == 1);
  assert(close(gradient.gradient_m[0], finite_difference, 1e-6));

  std::cout << "PASS: Quadra stats distributions, AR(1), MA(1), and AD gradient\n";
  return 0;
}
