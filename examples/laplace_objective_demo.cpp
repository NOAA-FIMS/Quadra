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
    RandomInterceptModel model({4.8, 5.1, 5.0, 4.9, 5.2});

    std::vector<double> theta = {4.7};
    std::vector<double> u0 = {0.0};

    auto result =
        quadra::evaluate_laplace_objective(
            model,
            theta,
            u0,
            model.parameters()
        );

    std::cout << "Laplace objective demo\n";
    std::cout << "converged = " << result.converged_m << "\n";
    std::cout << "logdet ok = " << result.logdet_ok_m << "\n";
    std::cout << "message = " << result.message_m << "\n";
    std::cout << "u_hat = " << result.u_hat_m[0] << "\n";
    std::cout << "joint objective = " << result.joint_objective_m << "\n";
    std::cout << "logdet(H_uu) = " << result.log_det_hessian_m << "\n";
    std::cout << "Laplace objective = " << result.laplace_objective_m << "\n";

    return 0;
}
