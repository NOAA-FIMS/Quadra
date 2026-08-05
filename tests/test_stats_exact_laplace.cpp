#include "../core/model/quadra_model.hpp"
#include "../include/quadra/stats.hpp"

#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

DECLARE_ADGRAPH();

class ThetaDependentRandomModel
    : public quadra::QuadraModel<ThetaDependentRandomModel> {
public:
  ThetaDependentRandomModel() {
    parameters_m.add("theta", 0.4, quadra::ParameterTransform::Identity,
                     false);
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
    const T u = parameters[1];
    return T(0.5) * theta * theta +
           T(0.5) * (T(1.0) + exp(theta)) * u * u;
  }

private:
  quadra::ParameterSet parameters_m;
};

class BranchingStructureModel
    : public quadra::QuadraModel<BranchingStructureModel> {
public:
  BranchingStructureModel() {
    parameters_m.add("theta", -1.0, quadra::ParameterTransform::Identity,
                     false);
    parameters_m.add("u0", 0.0, quadra::ParameterTransform::Identity, true);
    parameters_m.add("u1", 0.0, quadra::ParameterTransform::Identity, true);
  }
  std::vector<std::string> parameter_names_impl() const {
    return parameters_m.names();
  }
  const quadra::ParameterSet &parameters() const { return parameters_m; }
  template <typename T>
  T evaluate_impl(const std::vector<T> &p,
                  quadra::ModelReportContext &) const {
    T objective = T(0.5) * (p[1] * p[1] + p[2] * p[2]);
    if (quadra::value_of(p[0]) > 0.0) {
      objective += T(0.2) * p[1] * p[2];
    }
    return objective;
  }

private:
  quadra::ParameterSet parameters_m;
};

class TwoFixedEffectModel
    : public quadra::QuadraModel<TwoFixedEffectModel> {
public:
  TwoFixedEffectModel() {
    parameters_m.add("theta", 0.4, quadra::ParameterTransform::Identity,
                     false);
    parameters_m.add("scale", -0.1, quadra::ParameterTransform::Identity,
                     false);
    parameters_m.add("location", 0.2, quadra::ParameterTransform::Identity,
                     false);
    parameters_m.add("u", 0.0, quadra::ParameterTransform::Identity, true);
  }
  std::vector<std::string> parameter_names_impl() const {
    return parameters_m.names();
  }
  const quadra::ParameterSet &parameters() const { return parameters_m; }
  template <typename T>
  T evaluate_impl(const std::vector<T> &p,
                  quadra::ModelReportContext &) const {
    return T(0.5) * p[2] * p[2] +
           T(0.5) * (T(1.0) + exp(p[0]) + exp(p[1])) * p[3] * p[3];
  }

private:
  quadra::ParameterSet parameters_m;
};

int main() {
  ThetaDependentRandomModel model;
  quadra::LaplaceObjectiveOptions options;
  options.include_constant_m = true;
  options.newton_m.gradient_tolerance_m = 1e-12;

  quadra::stats::ExactLaplaceEvaluator<ThetaDependentRandomModel> evaluator(
      model, {0.4}, {0.0}, model.parameters(), options);
  const auto result = evaluator.evaluate({0.4}, {0.0});

  const double expected =
      0.4 + 0.5 * std::exp(0.4) / (1.0 + std::exp(0.4));
  if (!result.success || std::abs(result.gradient[0] - expected) > 1e-10) {
    std::cerr << "exact gradient mismatch: got " << result.gradient[0]
              << ", expected " << expected << "\n";
    return 1;
  }
  if (result.active_directions.size() != 1 ||
      result.active_directions[0] != 0) {
    std::cerr << "theta-dependent Hessian direction was not discovered\n";
    return 1;
  }
  if (result.objective.backend_m.backend !=
      quadra::laplace::LaplaceBackendKind::Diagonal) {
    std::cerr << "automatic diagonal backend was not selected\n";
    return 1;
  }

  const auto warm_result = evaluator.evaluate({0.4});
  if (!warm_result.success || evaluator.random_mode().size() != 1) {
    std::cerr << "stateful exact evaluator warm start failed\n";
    return 1;
  }

  quadra::stats::LaplaceEvaluator<ThetaDependentRandomModel> automatic(
      model, {0.0}, model.parameters(), options);
  const auto automatic_first = automatic.evaluate({0.4});
  const auto automatic_warm = automatic.evaluate({0.4});
  if (!automatic_first.converged_m || !automatic_first.logdet_ok_m ||
      automatic_first.backend_m.backend !=
          quadra::laplace::LaplaceBackendKind::Diagonal ||
      !automatic_first.structure_detected_m ||
      automatic_warm.structure_detected_m ||
      std::abs(automatic_first.laplace_objective_m - result.objective.laplace_objective_m) >
          1e-12) {
    std::cerr << "automatic persistent Laplace evaluator failed\n";
    return 1;
  }

  BranchingStructureModel branching_model;
  quadra::stats::LaplaceEvaluator<BranchingStructureModel> branching(
      branching_model, {0.0, 0.0}, branching_model.parameters(), options);
  const auto diagonal_branch = branching.evaluate({-1.0});
  const auto coupled_branch = branching.evaluate({1.0});
  if (diagonal_branch.backend_m.backend !=
          quadra::laplace::LaplaceBackendKind::Diagonal ||
      coupled_branch.backend_m.backend !=
          quadra::laplace::LaplaceBackendKind::Tridiagonal ||
      !coupled_branch.tape_rebuilt_m || branching.tape_rebuild_count() != 1) {
    std::cerr << "parameter-dependent tape topology was not rebuilt safely\n";
    return 1;
  }

  quadra::stats::ExactLaplaceEvaluator<BranchingStructureModel>
      branching_exact(branching_model, {-1.0}, {0.0, 0.0},
                      branching_model.parameters(), options);
  const auto rebuilt_exact = branching_exact.evaluate({1.0}, {0.0, 0.0});
  if (!rebuilt_exact.success ||
      branching_exact.hdot_tape_rebuild_count() != 1) {
    std::cerr << "persistent Hdot tape topology was not rebuilt safely\n";
    return 1;
  }

  quadra::laplace::ExactLaplaceGradientEngineOptions conservative_options;
  conservative_options.discover_active_directions = false;
  quadra::stats::ExactLaplaceEvaluator<BranchingStructureModel>
      conservative_exact(branching_model, {-1.0}, {0.0, 0.0},
                         branching_model.parameters(), options,
                         conservative_options);
  const auto conservative_result = conservative_exact.evaluate({-1.0});
  if (!conservative_result.success ||
      conservative_result.active_directions != std::vector<int>{0}) {
    std::cerr << "disabled direction discovery did not retain all directions\n";
    return 1;
  }

  TwoFixedEffectModel two_fixed_model;
  quadra::laplace::ExactLaplaceGradientEngineOptions worker_options;
  worker_options.hdot_workers = 2;
  quadra::stats::ExactLaplaceEvaluator<TwoFixedEffectModel> two_fixed(
      two_fixed_model, {0.4, -0.1, 0.2}, {0.0},
      two_fixed_model.parameters(), options, worker_options);
  const auto two_fixed_result = two_fixed.evaluate({0.4, -0.1, 0.2});
  if (!two_fixed_result.success ||
      two_fixed_result.active_directions != std::vector<int>({0, 1}) ||
      two_fixed.hdot_worker_count() != 2 ||
      two_fixed.hdot_shared_topology_owner_count() != 2 ||
      two_fixed.hdot_shared_operation_owner_count() != 2) {
    std::cerr << "multi-fixed Hdot workers or direction pruning failed\n";
    return 1;
  }
  worker_options.hdot_workers = 1;
  quadra::stats::ExactLaplaceEvaluator<TwoFixedEffectModel> serial_two_fixed(
      two_fixed_model, {0.4, -0.1, 0.2}, {0.0},
      two_fixed_model.parameters(), options, worker_options);
  const auto serial_two_fixed_result =
      serial_two_fixed.evaluate({0.4, -0.1, 0.2});
  if (!serial_two_fixed_result.success ||
      serial_two_fixed_result.gradient.size() !=
          two_fixed_result.gradient.size()) {
    std::cerr << "serial Hdot reference failed\n";
    return 1;
  }
  for (std::size_t i = 0; i < two_fixed_result.gradient.size(); ++i) {
    if (std::abs(two_fixed_result.gradient[i] -
                 serial_two_fixed_result.gradient[i]) > 1.0e-12) {
      std::cerr << "parallel Hdot gradient differs from serial reference\n";
      return 1;
    }
  }

  std::cout << "PASS: public stats exact Laplace gradient\n";
  return 0;
}
