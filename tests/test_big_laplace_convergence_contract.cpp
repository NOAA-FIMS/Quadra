#include "../examples/big/catch_at_age_laplace.cpp"

#include <cmath>
#include <iostream>
#include <vector>

int main() {
    using example::CatchAtAgeLaplaceModel;

    CatchAtAgeLaplaceModel model;

    quadra::ParameterSet params;

    // Fixed effects.
    params.add("log_R0", std::log(900.0), quadra::Transform::Identity, false);
    params.add("log_M", std::log(0.25), quadra::Transform::Identity, false);
    params.add("log_q", std::log(0.15), quadra::Transform::Identity, false);
    params.add("log_Fbar", std::log(0.18), quadra::Transform::Identity, false);
    params.add("sel50_raw", 0.0, quadra::Transform::Identity, false);
    params.add("log_sel_slope", std::log(1.25), quadra::Transform::Identity, false);
    params.add("log_sigma_index", std::log(0.20), quadra::Transform::Identity, false);
    params.add("log_sigma_catch", std::log(0.18), quadra::Transform::Identity, false);
    params.add("log_sigma_rec", std::log(0.35), quadra::Transform::Identity, false);

    // Random effects: 30 annual recruitment deviations.
    for (int y = 0; y < 30; ++y) {
        params.add("rec_dev_" + std::to_string(y + 1),
                   0.0,
                   quadra::Transform::Identity,
                   true);
    }

    auto opts = quadra::default_laplace_options();
    opts.gradient_tolerance = 1.0e-4;
    opts.max_iterations = 400;

    auto fit = quadra::optimize_lbfgs(model, params, opts);

    quadra::LaplaceEvaluator<CatchAtAgeLaplaceModel> evaluator(model, params, opts);
    const auto result = evaluator.evaluate(fit.par);

    std::vector<double> full = fit.par;
    full.insert(full.end(),
                result.random_effects_m.begin(),
                result.random_effects_m.end());

    const double direct_full = model(full);
    const double direct_diff = direct_full - result.joint_objective_m;

    bool ok = true;

    ok = ok && std::isfinite(fit.value);
    ok = ok && std::isfinite(result.laplace_objective_m);
    ok = ok && fit.grad_norm < 1.0e-4;
    ok = ok && result.random_gradient_norm_m < 1.0e-6;
    ok = ok && std::fabs(direct_diff) < 1.0e-8;

    if (!ok) {
        std::cerr << "FAIL: big Laplace black-box convergence contract failed\n";
        std::cerr << "  fit value: " << fit.value << "\n";
        std::cerr << "  Laplace objective: " << result.laplace_objective_m << "\n";
        std::cerr << "  fixed gradient norm: " << fit.grad_norm << "\n";
        std::cerr << "  random gradient norm: " << result.random_gradient_norm_m << "\n";
        std::cerr << "  direct_full_minus_reported: " << direct_diff << "\n";
        return 1;
    }

    std::cout << "PASS: big Laplace black-box convergence contract satisfied\n";
    std::cout << "  fit value: " << fit.value << "\n";
    std::cout << "  Laplace objective: " << result.laplace_objective_m << "\n";
    std::cout << "  fixed gradient norm: " << fit.grad_norm << "\n";
    std::cout << "  random gradient norm: " << result.random_gradient_norm_m << "\n";
    std::cout << "  direct_full_minus_reported: " << direct_diff << "\n";

    return 0;
}
