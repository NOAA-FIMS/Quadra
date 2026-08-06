#include "../../core/laplace/laplace_implicit_workspace.hpp"
#include "../../include/quadra/stats.hpp"
#include "catch_at_age_inference.hpp"
#include "catch_at_age_shared.hpp"

int main() {
  example::CatchAtAgeLaplaceModel model;
  quadra::ParameterSet parameters;
  quadra::ParameterVector optimizer_parameters;
  std::vector<double> fixed;
  fixed.reserve(10);

  auto add_fixed = [&](const std::string &name, const double value) {
    parameters.add(name, value, quadra::ParameterTransform::Identity, false);
    optimizer_parameters.add(
        {name, value, quadra::ParameterTransform::Identity, false});
    fixed.push_back(value);
  };

  add_fixed("log_R0", std::log(900.0));
  add_fixed("log_M", std::log(0.25));
  add_fixed("log_q", std::log(0.15));
  add_fixed("log_Fbar", std::log(0.18));
  add_fixed("sel50_raw", 0.0);
  add_fixed("log_sel_slope", std::log(1.25));
  add_fixed("log_sigma_index", std::log(0.20));
  add_fixed("log_sigma_catch", std::log(0.15));
  add_fixed("log_sigma_rec", std::log(0.35));
  add_fixed("log_comp_concentration", std::log(40.0));

  std::vector<double> random_initial(
      static_cast<std::size_t>(model.data.n_years), 0.0);
  for (int year = 0; year < model.data.n_years; ++year) {
    parameters.add("rec_dev_" + std::to_string(year + 1), 0.0,
                   quadra::ParameterTransform::Identity, true);
    optimizer_parameters.add({"rec_dev_" + std::to_string(year + 1), 0.0,
                              quadra::ParameterTransform::Identity, true});
  }

  quadra::LaplaceObjectiveOptions objective_options;
  objective_options.include_constant_m = true;
  objective_options.newton_m.gradient_tolerance_m = 1.0e-8;

  quadra::LaplaceOptions optimizer_options = quadra::default_laplace_options();
  optimizer_options.use_hutchinson_trace = false;
  optimizer_options.hessian_drop_tol = 0.0;
  const auto fit =
      quadra::optimize_lbfgs(model, optimizer_parameters, optimizer_options);

  if (!std::isfinite(fit.value) || fit.u_hat.size() != random_initial.size()) {
    std::cerr << "Big catch-at-age Laplace fit failed: " << fit.message << "\n";
    return 1;
  }

  quadra::stats::LaplaceEvaluator<example::CatchAtAgeLaplaceModel> evaluator(
      model, fit.u_hat, parameters, objective_options);
  const auto result = evaluator.evaluate(fit.par);
  const std::vector<double> &final_random_effects = result.u_hat_m;
  std::vector<std::string> fixed_parameter_names;
  fixed_parameter_names.reserve(fixed.size());
  const auto partition = quadra::partition_parameters(parameters);
  const auto parameter_names = parameters.names();
  for (const int index : partition.fixed_indices_m) {
    fixed_parameter_names.push_back(
        parameter_names[static_cast<std::size_t>(index)]);
  }

  quadra::laplace::ExactLaplaceGradientEngineOptions inference_engine_options;
  inference_engine_options.discover_active_directions = false;
  auto inference_objective_options = objective_options;
  // Covariance perturbations start extremely close to the fitted mode. Asking
  // Newton to improve beyond this scale can make its line search reject every
  // step because of floating-point noise.
  inference_objective_options.newton_m.gradient_tolerance_m = 1.0e-6;
  quadra::stats::ExactLaplaceEvaluator<example::CatchAtAgeLaplaceModel>
      covariance_evaluator(model, fit.par, final_random_effects, parameters,
                           inference_objective_options,
                           inference_engine_options);
  auto laplace_gradient_for_covariance =
      [&](const std::vector<double> &theta) -> std::vector<double> {
    const auto evaluation = covariance_evaluator.evaluate(theta);
    if (!evaluation.success) {
      throw std::runtime_error(
          "exact Laplace gradient evaluation failed during inference: " +
          evaluation.objective.message_m);
    }
    return evaluation.gradient;
  };

  const auto implicit_workspace = quadra::build_laplace_implicit_workspace(
      model, fit.par, final_random_effects, parameters);
  example::run_big_catch_at_age_inference(
      laplace_gradient_for_covariance, model, final_random_effects,
      fixed_parameter_names, fit.par, implicit_workspace);

  std::vector<double> full_parameters = fit.par;
  full_parameters.insert(full_parameters.end(), final_random_effects.begin(),
                         final_random_effects.end());
  const double direct_joint_objective = model(full_parameters);

  std::cout << "\nQuadra big fisheries example: catch-at-age Laplace model\n";
  std::cout << "years: " << model.data.n_years << "\n";
  std::cout << "ages: " << model.data.n_ages << "\n";
  std::cout << "fixed effects: " << fit.par.size() << "\n";
  std::cout << "random effects: " << final_random_effects.size() << "\n";
  std::cout << "converged: " << (fit.converged ? "yes" : "no") << "\n";
  std::cout << "message: " << fit.message << "\n";
  if (!fit.converged) {
    std::cout << "warning: this diagnostic example reached a finite optimizer "
                 "plateau; do not treat its estimates as a converged "
                 "assessment\n";
  }
  std::cout << "Laplace objective: " << fit.value << "\n";
  std::cout << "fixed gradient norm: " << fit.grad_norm << "\n";
  if (fit.fixed_gradient.size() == fit.fixed_gradient_names.size()) {
    std::cout << "fixed gradient components:\n";
    for (std::size_t i = 0; i < fit.fixed_gradient.size(); ++i) {
      std::cout << "  " << std::setw(24) << std::left
                << fit.fixed_gradient_names[i] << std::right << " "
                << fit.fixed_gradient[i] << "\n";
    }
  }
  std::cout << "joint objective: " << result.joint_objective_m << "\n";
  std::cout << "direct joint difference: "
            << direct_joint_objective - result.joint_objective_m << "\n";
  std::cout << "backend: "
            << quadra::laplace::ToString(result.backend_m.backend) << "\n";

  example::print_fixed_parameter_report(fit.par);
  example::print_objective_decomposition_report(model, fit.par, result,
                                                final_random_effects);
  print_actual_objective_path_decomposition(model, fit.par, result,
                                            final_random_effects);
  example::print_index_catch_diagnostics(model, fit.par, result,
                                         final_random_effects);
  example::print_age_composition_diagnostics(model, fit.par, result,
                                             final_random_effects);
  example::print_derived_quantity_report(model, fit.par, result,
                                         final_random_effects);

  // This is an executable diagnostics example, not a convergence contract.
  // A finite completed workflow is success; convergence status is reported
  // explicitly above and tested separately by optimizer contracts.
  return 0;
}
