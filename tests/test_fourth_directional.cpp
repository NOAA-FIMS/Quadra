#include "../core/autodiff.hpp"
#include "../core/laplace/laplace_exact_directional_curvature.hpp"
#include "../examples/big/catch_at_age_shared.hpp"

#include <cmath>
#include <iostream>
#include <vector>

struct CurvatureModel {
  quadra::ParameterSet parameters;
  CurvatureModel() {
    parameters.add("theta", 0.0, quadra::ParameterTransform::Identity, false);
    parameters.add("u", 0.0, quadra::ParameterTransform::Identity, true);
  }
  template <class Context> void initialize(Context &) {}
  template <class T> T evaluate(const std::vector<T> &p) const {
    const T residual = T(2.0) - p[0] - p[1];
    return T(2.5) * residual * residual + T(0.5) * p[1] * p[1];
  }
  template <class T, class Context>
  T evaluate(const std::vector<T> &p, Context &) const { return evaluate(p); }
};

struct ThetaCurvatureModel : CurvatureModel {
  template <class T> T evaluate(const std::vector<T> &p) const {
    using had::exp;
    return T(0.5) * p[0] * p[0] + T(0.5) * exp(p[0] * p[0]) * p[1] * p[1];
  }
  template <class T, class Context>
  T evaluate(const std::vector<T> &p, Context &) const { return evaluate(p); }
};

int main() {
  const std::vector<double> x{0.3, -0.2};
  const std::vector<double> d{1.5, -0.7};
  auto function = [](const auto &z) {
    using had::exp;
    using had::log;
    return exp(z[0] * z[1]) + log(z[0] + 2.0) + z[0] * z[0] * z[0] * z[0];
  };
  const auto got = had::evaluate_directional_derivatives4(function, x, d);
  auto univariate = [&](double t) {
    const double a = x[0] + t * d[0];
    const double b = x[1] + t * d[1];
    return std::exp(a * b) + std::log(a + 2.0) + std::pow(a, 4.0);
  };
  // Five-point central fourth derivative is only a loose independent check;
  // exact polynomial and exponential identities below provide tight checks.
  const double h = 2.0e-3;
  const double fd4 = (univariate(-2*h) - 4*univariate(-h) +
                      6*univariate(0.0) - 4*univariate(h) +
                      univariate(2*h)) / std::pow(h, 4.0);
  if (std::abs(got.fourth - fd4) > 2.0e-2) {
    std::cerr << "FAIL: fourth directional derivative mismatch: "
              << got.fourth << " vs " << fd4 << '\n';
    return 1;
  }

  example::CatchAtAgeLaplaceModel model(2, 3);
  std::vector<double> parameters{
      std::log(900.0), std::log(0.25), std::log(0.15), std::log(0.18), 0.0,
      std::log(1.25), std::log(0.20), std::log(0.15), std::log(0.35),
      std::log(40.0), 0.0, 0.0};
  std::vector<double> direction(parameters.size(), 0.0);
  direction[0] = 1.0;
  direction[10] = -0.25;
  const auto catch_result = had::evaluate_directional_derivatives4(
      [&model](const auto &p) { return model.evaluate(p); }, parameters,
      direction);
  if (!std::isfinite(catch_result.fourth)) {
    std::cerr << "FAIL: catch-at-age fourth derivative is not finite\n";
    return 1;
  }

  CurvatureModel curvature_model;
  const auto partition = quadra::partition_parameters(curvature_model.parameters);
  const auto curvature = quadra::laplace::exact_laplace_directional_curvature(
      curvature_model, std::vector<double>{0.0},
      std::vector<double>{5.0 / 3.0}, partition,
      Eigen::VectorXd::Ones(1));
  if (std::abs(curvature.curvature_m - 5.0 / 6.0) > 1.0e-11 ||
      std::abs(curvature.logdet_curvature_m) > 1.0e-11) {
    std::cerr << "FAIL: exact Laplace directional curvature mismatch\n";
    return 1;
  }
  ThetaCurvatureModel theta_curvature_model;
  const auto theta_curvature =
      quadra::laplace::exact_laplace_directional_curvature(
          theta_curvature_model, std::vector<double>{0.3},
          std::vector<double>{0.0},
          quadra::partition_parameters(theta_curvature_model.parameters),
          Eigen::VectorXd::Ones(1));
  if (std::abs(theta_curvature.joint_profile_curvature_m - 1.0) > 1.0e-11 ||
      std::abs(theta_curvature.logdet_curvature_m - 2.0) > 1.0e-10 ||
      std::abs(theta_curvature.curvature_m - 2.0) > 1.0e-10) {
    std::cerr << "FAIL: fourth-order logdet curvature mismatch\n";
    return 1;
  }
  std::cout << "PASS: fourth directional AD\n";
  return 0;
}
