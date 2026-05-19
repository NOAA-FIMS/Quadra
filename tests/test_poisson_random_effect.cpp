
#include "test_common.hpp"

DECLARE_ADGRAPH();

struct PoissonRandomIntercept {
    std::vector<double> y;

    template <typename T>
    T operator()(const std::vector<T>& p) const {
        T log_lambda = p[0];
        T log_sigma = p[1];
        T sigma = exp(log_sigma);

        T nll = 0.0;

        for (int i = 0; i < static_cast<int>(y.size()); ++i) {
            T u = p[2 + i];
            T eta = log_lambda + u;
            T lambda = exp(eta);

            // Poisson negative log likelihood, ignoring log(y!)
            nll += lambda - y[i] * eta;

            // Gaussian RE prior
            nll += 0.5 * (u / sigma) * (u / sigma) + log_sigma;
        }

        // weak penalties for fixed effects
        nll += 0.01 * log_lambda * log_lambda;
        nll += 0.01 * log_sigma * log_sigma;

        return nll;
    }
};

int main() {
    using namespace quadra;
    quadra_tests::print_banner("Testing PoissonRandomIntercept");

    std::vector<double> y = {3, 4, 2, 5, 7, 3, 4, 6};

    ParameterVector params;
    params.add({"log_lambda", std::log(4.0), ParameterTransform::Identity, false});
    params.add({"log_sigma", std::log(0.5), ParameterTransform::Identity, false});

    for (int i = 0; i < static_cast<int>(y.size()); ++i) {
        params.add({"u_" + std::to_string(i), 0.0, ParameterTransform::Identity, true});
    }

    PoissonRandomIntercept model{y};

    auto opts = quadra_tests::default_test_options();
    auto fit = optimize_lbfgs(model, params, opts);

    std::cout << "fit.value = " << fit.value << "\n";
    std::cout << "log_lambda_hat = " << fit.par[0] << "\n";
    std::cout << "log_sigma_hat  = " << fit.par[1] << "\n";

    if (!std::isfinite(fit.value)) {
        std::cerr << "FAIL: non-finite objective.\n";
        return 1;
    }

    std::cout << "PASS\n";
    return 0;
}
