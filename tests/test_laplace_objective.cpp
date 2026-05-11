#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "../core/laplace/laplace_objective.hpp"
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

        // Observation model: y_i ~ N(mu + u, 1), constants omitted.
        for (double yi : y_m) {
            Type r = Type(yi) - (mu + u);
            nll += Type(0.5) * r * r;
        }

        // Random-effect prior: u ~ N(0, 1), constants omitted.
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
    std::cout << "Testing Laplace objective\n";

    RandomInterceptModel model({4.8, 5.1, 5.0, 4.9, 5.2});

    std::vector<double> fixed = {4.7};
    std::vector<double> u0 = {0.0};

    quadra::LaplaceObjectiveOptions options;
    options.include_constant_m = true;
    options.newton_m.gradient_tolerance_m = 1e-10;
    options.newton_m.step_tolerance_m = 1e-12;

    auto result =
        quadra::evaluate_laplace_objective(
            model,
            fixed,
            u0,
            model.parameters(),
            options
        );

    // Analytic solution:
    // u_hat = (sum(y) - n*mu) / (n + 1)
    //       = (25 - 5*4.7) / 6
    //       = 0.25.
    //
    // Joint objective at u_hat:
    // residuals are y_i - (4.7 + 0.25) = y_i - 4.95
    // residuals: -0.15, 0.15, 0.05, -0.05, 0.25
    // sum residual^2 = 0.1125
    // obs contribution = 0.05625
    // prior contribution = 0.5 * 0.25^2 = 0.03125
    // joint = 0.0875.
    //
    // H_uu = n + 1 = 6.
    //
    // Laplace objective with constant:
    // joint + 0.5 * log(6) - 0.5 * log(2*pi).
    const double expected_u_hat = 0.25;
    const double expected_joint = 0.0875;
    const double expected_logdet = std::log(6.0);
    const double expected_laplace =
        expected_joint + 0.5 * expected_logdet - 0.5 * std::log(2.0 * M_PI);

    std::cout << "converged = " << result.converged_m << "\n";
    std::cout << "message = " << result.message_m << "\n";
    std::cout << "u_hat = " << result.u_hat_m[0] << "\n";
    std::cout << "joint objective = " << result.joint_objective_m << "\n";
    std::cout << "logdet = " << result.log_det_hessian_m << "\n";
    std::cout << "laplace objective = " << result.laplace_objective_m << "\n";
    std::cout << "gradient norm random = " << result.gradient_norm_random_m << "\n";

    if (!result.converged_m) {
        std::cerr << "FAIL: Newton optimization did not converge\n";
        return 1;
    }

    if (!result.logdet_ok_m) {
        std::cerr << "FAIL: logdet computation failed\n";
        return 1;
    }

    if (std::abs(result.u_hat_m[0] - expected_u_hat) > 1e-10) {
        std::cerr << "FAIL: u_hat mismatch\n";
        return 1;
    }

    if (std::abs(result.joint_objective_m - expected_joint) > 1e-10) {
        std::cerr << "FAIL: joint objective mismatch\n";
        return 1;
    }

    if (std::abs(result.log_det_hessian_m - expected_logdet) > 1e-10) {
        std::cerr << "FAIL: logdet mismatch\n";
        return 1;
    }

    if (std::abs(result.laplace_objective_m - expected_laplace) > 1e-10) {
        std::cerr << "FAIL: Laplace objective mismatch\n";
        return 1;
    }

    if (result.gradient_norm_random_m > 1e-10) {
        std::cerr << "FAIL: random gradient norm too large\n";
        return 1;
    }

    std::cout << "PASS\n";
    return 0;
}
