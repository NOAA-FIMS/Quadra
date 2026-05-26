#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "../core/laplace/laplace_profile.hpp"
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
    std::cout << "Testing Laplace profiling wrappers\n";

    RandomInterceptModel model({4.8, 5.1, 5.0, 4.9, 5.2});

    auto component_profile =
        quadra::profile_laplace_components(
            model,
            {4.7},
            {0.0},
            model.parameters()
        );

    auto optimizer_profile =
        quadra::profile_exact_lbfgs(
            model,
            {4.0},
            {0.0},
            model.parameters()
        );

    std::cout << "component total ms = " << component_profile.total_ms_m << "\n";
    std::cout << "newton ms = " << component_profile.random_newton_ms_m << "\n";
    std::cout << "hessian ms = " << component_profile.random_hessian_ms_m << "\n";
    std::cout << "laplace objective ms = " << component_profile.laplace_objective_ms_m << "\n";
    std::cout << "exact gradient ms = " << component_profile.exact_gradient_ms_m << "\n";

    std::cout << "exact lbfgs ms = " << optimizer_profile.exact_lbfgs_ms_m << "\n";
    std::cout << "exact lbfgs iter = " << optimizer_profile.exact_lbfgs_iterations_m << "\n";

    if (!(component_profile.total_ms_m >= 0.0)) return 1;
    if (!(component_profile.random_newton_ms_m >= 0.0)) return 1;
    if (!(component_profile.random_hessian_ms_m >= 0.0)) return 1;
    if (!(component_profile.exact_gradient_ms_m >= 0.0)) return 1;
    if (!optimizer_profile.exact_lbfgs_converged_m) return 1;

    std::cout << "PASS\n";
    return 0;
}
