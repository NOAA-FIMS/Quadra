#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "../core/laplace/laplace_exact_objective_gradient.hpp"
#include "../core/model/quadra_model.hpp"

DECLARE_ADGRAPH();

class RandomInterceptModel : public quadra::QuadraModel<RandomInterceptModel> {
public:
    explicit RandomInterceptModel(std::vector<double> y)
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
    Type evaluate_impl(const std::vector<Type>& p, quadra::ModelReportContext&) const {
        Type mu = p[0];
        Type u = p[1];
        Type nll = Type(0.0);

        for (double yi : y_m) {
            Type r = Type(yi) - (mu + u);
            nll += Type(0.5) * r * r;
        }

        nll += Type(0.5) * u * u;
        return nll;
    }

private:
    std::vector<double> y_m;
    quadra::ParameterSet parameters_m;
};

int main() {
    RandomInterceptModel model({4.8, 5.1, 5.0, 4.9, 5.2});

    auto result =
        quadra::evaluate_laplace_exact_objective_gradient(
            model,
            {4.7},
            {0.0},
            model.parameters()
        );

    std::cout << "gradient_mu = " << result.gradient_fixed_m[0] << "\n";
    std::cout << "objective = " << result.laplace_objective_m << "\n";

    if (!result.converged_m || !result.logdet_ok_m) return 1;
    if (std::abs(result.gradient_fixed_m[0] - (-0.25)) > 1e-12) return 1;

    std::cout << "PASS\n";
    return 0;
}
