#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "../core/laplace/laplace_exact_gradient.hpp"
#include "../core/model/quadra_model.hpp"
#include "../core/model/parameter_partition.hpp"

DECLARE_ADGRAPH();

class RandomInterceptModel : public quadra::QuadraModel<RandomInterceptModel> {
public:
    RandomInterceptModel(std::vector<double> y)
        : y_m(std::move(y))
    {
        parameters_m.add("mu", 0.0, quadra::ParameterTransform::Identity, false);
        parameters_m.add("u", 0.0, quadra::ParameterTransform::Identity, true);
    }

    std::vector<std::string> parameter_names_impl() const {
        return parameters_m.names();
    }

    const quadra::ParameterSet& parameters() const {
        return parameters_m;
    }

    template <typename Type>
    Type evaluate_impl(
        const std::vector<Type>& p,
        quadra::ModelReportContext& ctx
    ) const {
        Type mu = p[0];
        Type u = p[1];

        Type nll = Type(0.0);

        for (double yi : y_m) {
            Type r = Type(yi) - (mu + u);
            nll += Type(0.5) * r * r;
        }

        nll += Type(0.5) * u * u;

        ctx.report("mu", mu);
        ctx.report("u_hat", u);
        ctx.adreport("mean_prediction", mu + u);

        return nll;
    }

private:
    std::vector<double> y_m;
    quadra::ParameterSet parameters_m;
};

int main() {
    std::cout << "Testing exact Laplace envelope gradient\n";

    RandomInterceptModel model({4.8, 5.1, 5.0, 4.9, 5.2});

    std::vector<double> fixed = {4.7};
    std::vector<double> u0 = {0.0};

    quadra::LaplaceExactGradientOptions exact_options;
    exact_options.objective_m.include_constant_m = true;
    exact_options.objective_m.newton_m.gradient_tolerance_m = 1e-12;
    exact_options.objective_m.newton_m.step_tolerance_m = 1e-14;

    auto exact =
        quadra::evaluate_laplace_exact_gradient(
            model,
            fixed,
            u0,
            model.parameters(),
            exact_options
        );

    quadra::LaplaceFixedGradientOptions fd_options;
    fd_options.relative_step_m = 1e-6;
    fd_options.absolute_step_m = 1e-7;
    fd_options.use_central_difference_m = true;
    fd_options.objective_m = exact_options.objective_m;

    auto fd =
        quadra::evaluate_laplace_fixed_gradient(
            model,
            fixed,
            u0,
            model.parameters(),
            fd_options
        );

    // For this model H_uu = n + 1 is constant wrt mu, so
    // d/dmu logdet(H_uu) = 0. The envelope gradient is therefore the full
    // exact Laplace gradient.
    //
    // At mu = 4.7, u_hat = 0.25:
    // d/dmu f(mu, u_hat) = sum_i (mu + u_hat - y_i)
    //                    = 5 * 4.95 - 25
    //                    = -0.25.
    const double expected_gradient = -0.25;
    const double expected_u_hat = 0.25;

    std::cout << "converged = " << exact.converged_m << "\n";
    std::cout << "logdet ok = " << exact.logdet_ok_m << "\n";
    std::cout << "u_hat = " << exact.u_hat_m[0] << "\n";
    std::cout << "exact gradient_mu = " << exact.gradient_fixed_m[0] << "\n";
    std::cout << "finite-difference gradient_mu = " << fd.gradient_fixed_m[0] << "\n";
    std::cout << "exact gradient norm = " << exact.gradient_norm_m << "\n";

    if (!exact.converged_m || !exact.logdet_ok_m) {
        std::cerr << "FAIL: exact gradient base objective failed\n";
        return 1;
    }

    if (std::abs(exact.u_hat_m[0] - expected_u_hat) > 1e-10) {
        std::cerr << "FAIL: u_hat mismatch\n";
        return 1;
    }

    if (std::abs(exact.gradient_fixed_m[0] - expected_gradient) > 1e-12) {
        std::cerr << "FAIL: exact gradient mismatch\n";
        return 1;
    }

    if (std::abs(exact.gradient_fixed_m[0] - fd.gradient_fixed_m[0]) > 1e-5) {
        std::cerr << "FAIL: exact gradient does not match finite difference\n";
        return 1;
    }

    std::cout << "PASS\n";
    return 0;
}
