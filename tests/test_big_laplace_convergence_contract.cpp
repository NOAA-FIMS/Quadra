#define QUADRA_BIG_LAPLACE_NO_MAIN
#include "../examples/big/catch_at_age_laplace.cpp"

#include <cmath>
#include <iostream>

namespace example
{

    inline quadra::ParameterVector make_big_laplace_parameter_vector()

    {

        quadra::ParameterVector param_vector;

        param_vector.add(quadra::Parameter("log_R0", std::log(900.0), quadra::ParameterTransform::Identity, false));

        param_vector.add(quadra::Parameter("log_M", std::log(0.25), quadra::ParameterTransform::Identity, false));

        param_vector.add(quadra::Parameter("log_q", std::log(0.15), quadra::ParameterTransform::Identity, false));

        param_vector.add(quadra::Parameter("log_Fbar", std::log(0.18), quadra::ParameterTransform::Identity, false));

        param_vector.add(quadra::Parameter("sel50_raw", 0.0, quadra::ParameterTransform::Identity, false));

        param_vector.add(quadra::Parameter("log_sel_slope", std::log(1.25), quadra::ParameterTransform::Identity, false));

        param_vector.add(quadra::Parameter("log_sigma_index", std::log(0.20), quadra::ParameterTransform::Identity, false));

        param_vector.add(quadra::Parameter("log_sigma_catch", std::log(0.18), quadra::ParameterTransform::Identity, false));

        param_vector.add(quadra::Parameter("log_sigma_rec", std::log(0.35), quadra::ParameterTransform::Identity, false));

        for (int y = 0; y < 30; ++y)
        {

            param_vector.add(quadra::Parameter("rec_dev_" + std::to_string(y + 1),

                                               0.0,

                                               quadra::ParameterTransform::Identity,

                                               true));
        }

        return param_vector;
    }

} // namespace example

int main()
{
    example::CatchAtAgeLaplaceModel model;
    auto params = example::make_big_laplace_parameter_vector();

    quadra::LaplaceOptions opts = quadra::default_laplace_options();

    auto fit = quadra::optimize_lbfgs(model, params, opts);

    const bool ok =
        std::isfinite(fit.value) &&
        fit.value < 0.0 &&
        fit.value > -500.0 &&
        fit.grad_norm < 1.0e-3;

    if (!ok)
    {
        std::cerr << "FAIL: big Laplace black-box convergence contract failed\n";
        std::cerr << "  fit value: " << fit.value << "\n";
        std::cerr << "  fixed gradient norm: " << fit.grad_norm << "\n";
        return 1;
    }

    std::cout << "PASS: big Laplace black-box convergence contract satisfied\n";
    std::cout << "  fit value: " << fit.value << "\n";
    std::cout << "  fixed gradient norm: " << fit.grad_norm << "\n";
    return 0;
}
