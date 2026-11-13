#include "../core/autodiff/model_gradient.hpp"
#include "../core/model/quadra_model.hpp"
#include "../core/optimizer.hpp"
#include <iostream>
DECLARE_ADGRAPH();
struct Quadratic : quadra::QuadraModel<Quadratic> {
  std::vector<std::string> parameter_names_impl() const { return {"a", "b"}; }
  template <class T> T operator()(const std::vector<T> &x) const {
    const T a = x[0] - T(2), b = x[1] + T(1);
    return T(0.5) * a * a + T(2) * b * b;
  }
  template <class T>
  T evaluate_impl(const std::vector<T> &x, quadra::ModelReportContext &) const {
    return (*this)(x);
  }
};
struct MustNotEvaluate {
  template <class T> T operator()(const std::vector<T> &) const {
    throw std::runtime_error("Unexpected inner mode evaluation");
  }
};
int main() {
  try {
    quadra::ParameterVector p;
    p.add(quadra::Parameter("a", -3));
    p.add(quadra::Parameter("b", 4));
    had::ADGraph graph;
    MustNotEvaluate none;
    const auto mode = quadra::solve_random_effects_laplace(
        none, p, Eigen::VectorXd(2), {0, 1}, {}, graph);
    if (!mode.empty())
      throw std::runtime_error("Nonempty latent mode");
    Quadratic model;
    const auto gradient = quadra::evaluate_gradient(model, {-3, 4});
    if (std::abs(gradient.objective_value_m - 62.5) > 1e-12 ||
        std::abs(gradient.gradient_m[0] + 5) > 1e-12 ||
        std::abs(gradient.gradient_m[1] - 20) > 1e-12)
      throw std::runtime_error("First-order gradient incorrect");
    const auto fit = quadra::optimize_lbfgs(
        model, p, quadra::default_laplace_options(), 50, 1e-7);
    if (!fit.converged || std::abs(fit.par.at(0) - 2) > 1e-6 ||
        std::abs(fit.par.at(1) + 1) > 1e-6 ||
        std::abs(fit.laplace_logdet) > 1e-12 ||
        std::abs(fit.laplace_constant) > 1e-12 ||
        std::abs(fit.value - fit.joint_objective) > 1e-12)
      throw std::runtime_error("All-fixed public optimization incorrect");
    std::cout
        << "Fixed-only first-order evaluation and public optimizer: PASS\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
