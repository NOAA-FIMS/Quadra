#include "catch_at_age_shared.hpp"
#include "catch_at_age_inference.hpp"
#include "../../core/laplace/laplace_implicit_workspace.hpp"

int main()
{
    example::CatchAtAgeLaplaceModel model;

    quadra::ParameterSet parameters;
    quadra::ParameterVector param_vector;

    std::vector<double> fixed;
    fixed.reserve(9);

    auto add_fixed = [&](const std::string &name, double value)
    {
        parameters.add(name, value);
        param_vector.add({name, value, quadra::ParameterTransform::Identity, false});
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

    std::vector<double> random_initial;
    random_initial.reserve(static_cast<std::size_t>(model.data.n_years));

    for (int y = 0; y < model.data.n_years; ++y)
    {
        const std::string rec_name = "rec_dev_" + std::to_string(y + 1);

        parameters.add(
            rec_name,
            0.0,
            quadra::ParameterTransform::Identity,
            true);

        param_vector.add({rec_name,
                          0.0,
                          quadra::ParameterTransform::Identity,
                          true});

        random_initial.push_back(0.0);
    }

    quadra::LaplaceEvaluator<example::CatchAtAgeLaplaceModel> evaluator(
        model,
        parameters,
        random_initial);

    quadra::LaplaceOptions opts;
    opts.use_hutchinson_trace = false; // small example: exact deterministic trace
    opts.hessian_drop_tol = 0.0;
    
    // quadra_big_example_optimize_lbfgs_enabled
    auto fit = optimize_lbfgs(model, param_vector, opts);

    // quadra_big_example_final_eval_uses_fit_par
    const auto result = evaluator.evaluate(fit.par);

    const std::vector<int> final_fixed_idx = build_fixed_index(param_vector);
    const std::vector<int> final_random_idx = build_random_index(param_vector);
    had::ADGraph final_mode_graph;
    Eigen::VectorXd final_x(fit.par.size());
    for (int i = 0; i < final_x.size(); ++i)
        final_x[i] = fit.par[static_cast<size_t>(i)];

    const std::vector<double> final_random_effects =
        solve_random_effects_laplace(
            model,
            param_vector,
            final_x,
            final_fixed_idx,
            final_random_idx,
            final_mode_graph);


    // ========================================
    // Quadra inference integration
    // ========================================
    std::vector<std::string> fixed_parameter_names;
    fixed_parameter_names.reserve(static_cast<std::size_t>(final_x.size()));

    for (const int idx : final_fixed_idx)
    {
        fixed_parameter_names.push_back(param_vector.params[static_cast<std::size_t>(idx)].name);
    }

    auto laplace_objective_for_covariance =
        [&](const std::vector<double>& theta) -> double
        {
            quadra::LaplaceEvaluator<example::CatchAtAgeLaplaceModel> cov_evaluator(
                model,
                parameters,
                random_initial);

            const auto cov_result = cov_evaluator.evaluate(theta);
            return cov_result.laplace_objective_m;
        };

    // v1 derived placeholders: these scalar functions are constant wrappers
    // around final reported values. A later patch should replace these with
    // theta-dependent derived quantity functions for meaningful SEs.
    const auto implicit_workspace =
        quadra::build_laplace_implicit_workspace(
            model,
            fit.par,
            random_initial,
            parameters);

    example::run_big_catch_at_age_inference(
        laplace_objective_for_covariance,
        model,
        final_random_effects,
        fixed_parameter_names,
        fit.par,
        implicit_workspace);


    std::vector<double> full_par = fit.par;
    full_par.insert(full_par.end(),
                    final_random_effects.begin(),
                    final_random_effects.end());

    const double direct_full_objective = model(full_par);

    std::cout << "\nFull direct model objective check\n";
    std::cout << "fixed parameter count: " << fit.par.size() << "\n";
    std::cout << "random effect count: " << final_random_effects.size() << "\n";
    std::cout << "full parameter count: " << full_par.size() << "\n";
    std::cout << "model(fixed + random): " << direct_full_objective << "\n";
    std::cout << "result.joint_objective_m: " << result.joint_objective_m << "\n";
    std::cout << "direct_full_minus_reported: "
              << direct_full_objective - result.joint_objective_m << "\n";

    const double value = result.laplace_objective_m;

    std::cout << "Quadra big fisheries example: catch-at-age Laplace model\n";
    std::cout << "Final summary evaluated at optimize_lbfgs fit.par\n";
    std::cout << "years: " << model.data.n_years << "\n";
    std::cout << "ages: " << model.data.n_ages << "\n";
    std::cout << "fixed effects: " << fixed.size() << "\n";
    std::cout << "random effects: " << random_initial.size() << "\n";
    std::cout << "optimized fixed effects: " << fit.par.size() << "\n";

    std::cout << "\n";
    std::cout << "Optimizer summary\n";
    std::cout << "fit value: " << fit.value << "\n";
    std::cout << "fit iterations: " << fit.iterations << "\n";
    std::cout << "optimizer-reported fixed gradient norm: " << fit.grad_norm << "\n";
    example::print_fixed_parameter_report(fit.par);

    std::cout << "\n";
    std::cout << "Final Laplace re-evaluation\n";
    std::cout << "joint objective: " << result.joint_objective_m << "\n";
    std::cout << "log det Hessian: " << result.log_det_hessian_m << "\n";
    std::cout << "random gradient norm: " << result.gradient_norm_random_m << "\n";
    std::cout << "Laplace objective: " << result.laplace_objective_m << "\n";
    example::print_objective_decomposition_report(model, fit.par, result, final_random_effects);
    print_actual_objective_path_decomposition(model, fit.par, result, final_random_effects);
    example::print_index_catch_diagnostics(model, fit.par, result, final_random_effects);
    example::print_age_composition_diagnostics(model, fit.par, result, final_random_effects);
    example::print_derived_quantity_report(model, fit.par, result, final_random_effects);
    return 0;
}


