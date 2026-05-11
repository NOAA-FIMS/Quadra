#include <iostream>
#include <vector>

#include "../core/laplace/laplace_implicit_derivatives.hpp"
#include "../core/model/quadra_model.hpp"

DECLARE_ADGRAPH();

class DemoModel :
    public quadra::QuadraModel<DemoModel> {
public:
    DemoModel(std::vector<double> y)
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
    DemoModel model({4.8, 5.1, 5.0, 4.9, 5.2});

    auto result =
        quadra::evaluate_laplace_implicit_derivatives(
            model,
            {4.7},
            {0.0},
            model.parameters()
        );

    std::cout << "du/dtheta:\n";
    std::cout << result.du_dtheta_m << "\n";

    return 0;
}
