#include "core/autodiff/model_gradient.hpp"
#include "core/optimizer.hpp"
#include "swordfish/optimizer_scaling.hpp"
#include <chrono>
#include <filesystem>
#include <iostream>
DECLARE_ADGRAPH();

/// Analytic correlated quadratic used to check optimizer scaling independently
/// of fish data.
struct Quadratic : quadra::QuadraModel<Quadratic> {
  /// Stable coordinate labels for matrix identity checks.
  std::vector<std::string> parameter_names_impl() const { return {"a", "b"}; }
  /// Positive quadratic with minimum at zero and known correlated curvature.
  template <class T> T operator()(const std::vector<T> &x) const {
    return T(.5) *
           (T(4) * x[0] * x[0] + T(2) * x[0] * x[1] + T(2) * x[1] * x[1]);
  }
  /// Expose the analytic target through the standard model interface.
  template <class T>
  T evaluate_impl(const std::vector<T> &x, quadra::ModelReportContext &) const {
    return (*this)(x);
  }
};

/// Check objective/gradient chain rule, optimum recovery and metric label/SPD
/// rejection.
int main() {
  try {
    Quadratic model;
    Eigen::Matrix2d factor;
    factor << 2, 0, .3, .5;
    swordfish::OptimizerScaling<Quadratic> scaled(model, {.3, -.2}, factor);
    const std::vector<double> z{.1, -.1};
    const auto x = scaled.expand(z);
    const auto original = quadra::evaluate_gradient(model, x),
               gradient = quadra::evaluate_gradient(scaled, z);
    if (std::abs(original.objective_value_m - gradient.objective_value_m) >
        1e-12)
      throw std::runtime_error("Scaling changed objective");
    for (int j = 0; j < 2; ++j) {
      double expected = 0;
      for (int i = 0; i < 2; ++i)
        expected += factor(i, j) * original.gradient_m[i];
      if (std::abs(expected - gradient.gradient_m[j]) > 1e-12)
        throw std::runtime_error("Scaling chain rule");
    }
    quadra::ParameterVector parameters;
    parameters.add(quadra::Parameter("a", 0.));
    parameters.add(quadra::Parameter("b", 0.));
    const auto fit = quadra::optimize_lbfgs(
        scaled, parameters, quadra::default_laplace_options(), 100, 1e-9);
    const auto recovered = scaled.expand(fit.par);
    if (!fit.converged || std::hypot(recovered[0], recovered[1]) > 1e-8)
      throw std::runtime_error("Scaling changed optimum");
    const auto path =
        std::filesystem::temp_directory_path() /
        ("quadra-metric-" +
         std::to_string(
             std::chrono::steady_clock::now().time_since_epoch().count()) +
         ".csv");
    {
      std::ofstream out(path);
      out << "quantity,a,b\na,4,1\nb,1,2\n";
    }
    const auto matrix =
        swordfish::read_optimizer_metric(path.string(), {"a", "b"});
    if (matrix(0, 1) != 1)
      throw std::runtime_error("Metric read value");
    for (int mode = 0; mode < 3; ++mode) {
      {
        std::ofstream out(path);
        if (mode == 0)
          out << "quantity,b,a\na,4,1\nb,1,2\n";
        if (mode == 1)
          out << "quantity,a,b\na,1,2\nb,2,1\n";
        if (mode == 2)
          out << "quantity,a,b\na,4,1\n";
      }
      bool rejected = false;
      try {
        (void)swordfish::read_optimizer_metric(path.string(), {"a", "b"});
      } catch (const std::exception &) {
        rejected = true;
      }
      if (!rejected)
        throw std::runtime_error("Invalid metric accepted");
    }
    std::filesystem::remove(path);
    std::cout << "Optimizer scaling target, full chain rule, optimum and "
                 "labeled SPD metric: PASS\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
