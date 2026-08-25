#include "../core/model/quadra_model.hpp"
#include "../include/quadra/stats.hpp"

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

DECLARE_ADGRAPH();

class CurvatureModel : public quadra::QuadraModel<CurvatureModel> {
public:
  CurvatureModel() {
    parameters_m.add("theta", 0.4, quadra::ParameterTransform::Identity, false);
    parameters_m.add("u", 0.0, quadra::ParameterTransform::Identity, true);
  }

  std::vector<std::string> parameter_names_impl() const {
    return parameters_m.names();
  }
  const quadra::ParameterSet &parameters() const { return parameters_m; }

  template <typename T>
  T evaluate_impl(const std::vector<T> &parameters,
                  quadra::ModelReportContext &) const {
    const T theta = parameters[0];
    const T random = parameters[1];
    return T(0.5) * theta * theta +
           T(0.5) * (T(1.0) + exp(theta)) * random * random;
  }

private:
  quadra::ParameterSet parameters_m;
};

int main() {
  CurvatureModel model;
  quadra::LaplaceObjectiveOptions objective_options;
  objective_options.include_constant_m = true;
  objective_options.newton_m.gradient_tolerance_m = 1.0e-12;

  quadra::stats::ExactLaplaceEvaluator<CurvatureModel> evaluator(
      model, {0.4}, {0.0}, model.parameters(), objective_options);
  const double initial_objective =
      evaluator.evaluate({0.4}).objective.laplace_objective_m;

  quadra::stats::LaplaceOptimizerOptions optimizer_options;
  optimizer_options.gradient_tolerance = 1.0e-10;
  optimizer_options.step_tolerance = 1.0e-13;
  const auto result =
      quadra::stats::optimize_laplace(evaluator, {0.4}, optimizer_options);

  if (!result.converged || result.gradient.size() != 1 ||
      std::abs(result.gradient[0]) > 1.0e-9 ||
      !(result.objective < initial_objective) || result.history.empty() ||
      result.objective_evaluations == 0 || result.exact_evaluations == 0 ||
      evaluator.objective_tape_rebuild_count() != 0 ||
      evaluator.hdot_tape_rebuild_count() != 0) {
    std::cerr << "stateful exact Laplace optimizer failed: " << result.message
              << " gradient="
              << (result.gradient.empty() ? std::nan("") : result.gradient[0])
              << "\n";
    return 1;
  }

  const double theta = result.fixed[0];
  const double expected_gradient =
      theta + 0.5 * std::exp(theta) / (1.0 + std::exp(theta));
  if (std::abs(expected_gradient) > 1.0e-9 ||
      std::abs(result.random_mode[0]) > 1.0e-12) {
    std::cerr << "optimizer did not include the logdet derivative\n";
    return 1;
  }

  const auto convenience =
      quadra::stats::optimize_laplace(model, {0.4}, {0.0}, model.parameters(),
                                      optimizer_options, objective_options);
  if (!convenience.converged ||
      std::abs(convenience.fixed[0] - result.fixed[0]) > 1.0e-9) {
    std::cerr << "automatic optimizer construction failed\n";
    return 1;
  }

  quadra::stats::ExactLaplaceEvaluator<CurvatureModel> budget_evaluator(
      model, {0.4}, {0.0}, model.parameters(), objective_options);
  auto budget_options = optimizer_options;
  budget_options.max_evaluations = 1;
  const auto budget_result =
      quadra::stats::optimize_laplace(budget_evaluator, {0.4}, budget_options);
  if (budget_result.converged ||
      budget_result.message.find("budget exhausted") == std::string::npos) {
    std::cerr << "evaluation budget was not enforced\n";
    return 1;
  }

  quadra::stats::ExactLaplaceEvaluator<CurvatureModel> hybrid_evaluator(
      model, {0.4}, {0.0}, model.parameters(), objective_options);
  const auto hybrid = quadra::stats::optimize_laplace_hybrid(
      hybrid_evaluator, {0.4}, 0, 0.3, optimizer_options);
  if (!hybrid.converged || hybrid.approximate_evaluations == 0 ||
      hybrid.exact_evaluations == 0 || hybrid.exact_switch_iteration < 0 ||
      std::abs(hybrid.gradient[0]) > 1.0e-9) {
    std::cerr << "hybrid exact-switch convergence contract failed: "
              << hybrid.message << '\n';
    return 1;
  }

  std::cout << "PASS: public exact and hybrid Laplace optimizers\n";
  return 0;
}
