#include <cmath>
#include <iostream>
#include <vector>

#include "../core/laplace/laplace_implicit_derivatives.hpp"
#include "../core/model/quadra_model.hpp"

DECLARE_ADGRAPH();

class RandomInterceptModel :
    public quadra::QuadraModel<RandomInterceptModel> {
public:
    RandomInterceptModel(std::vector<double> y)
        : y_m(std::move(y))
    {
        parameters_m.add(
            "mu",
            0.0,
            quadra::ParameterTransform::Identity,
            false
        );

        parameters_m.add(
            "u",
            0.0,
            quadra::ParameterTransform::Identity,
            true
        );
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
        quadra::ModelReportContext&
    ) const {
        Type mu = p[0];
        Type u = p[1];

        Type nll = 0.0;

        for (double yi : y_m) {
            Type r = yi - (mu + u);
            nll += 0.5 * r * r;
        }

        nll += 0.5 * u * u;

        return nll;
    }

private:
    std::vector<double> y_m;
    quadra::ParameterSet parameters_m;
};

int main() {
    std::cout << "Testing implicit derivatives\n";

    RandomInterceptModel model({4.8, 5.1, 5.0, 4.9, 5.2});

    std::vector<double> theta = {4.7};
    std::vector<double> u0 = {0.0};

    auto result =
        quadra::evaluate_laplace_implicit_derivatives(
            model,
            theta,
            u0,
            model.parameters()
        );

    std::cout << "success = " << result.success_m << "\n";
    std::cout << "message = " << result.message_m << "\n";

    std::cout << "H_u_theta:\n";
    std::cout << result.H_u_theta_m << "\n";

    std::cout << "du_dtheta:\n";
    std::cout << result.du_dtheta_m << "\n";

    // Analytic:
    //
    // gradient_u = n(mu + u) - sum(y) + u
    //
    // d/du gradient_u = n + 1 = 6
    // d/dmu gradient_u = n = 5
    //
    // du/dmu = -5/6
    const double expected = -5.0 / 6.0;

    if (!result.success_m) {
        std::cerr << "FAIL: implicit derivative solve failed\n";
        return 1;
    }

    if (std::abs(result.du_dtheta_m(0,0) - expected) > 1e-4) {
        std::cerr << "FAIL: implicit derivative mismatch\n";
        return 1;
    }

    std::cout << "PASS\n";
    return 0;
}
