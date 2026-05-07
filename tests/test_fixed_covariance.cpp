#include "test_common.hpp"
#include "../core/inference/covariance.hpp"

struct CurvatureDependsOnTheta {
    int n_random;

    template <typename T>
    T operator()(const std::vector<T>& p) const {
        T theta = p[0];
        T nll = 0.5 * theta * theta;

        for (int i = 0; i < n_random; ++i) {
            T u = p[1 + i];
            nll += 0.5 * exp(theta) * u * u;
        }

        return nll;
    }
};

int main() {
    using namespace quadra;

    quadra_tests::print_banner("Testing fixed-effect covariance");

    int n_random = 10;

    ParameterVector params;
    params.add({"theta", 0.0, Transform::Identity, false});

    for (int i = 0; i < n_random; ++i) {
        params.add({"u_" + std::to_string(i), 0.0, Transform::Identity, true});
    }

    CurvatureDependsOnTheta model{n_random};

    LaplaceOptions lopts = quadra_tests::default_test_options();
    lopts.use_hutchinson_trace = false;

    auto fit = optimize_lbfgs(model, params, lopts);

    CovarianceOptions copts;
    copts.fd_step = 1e-4;
    copts.verbose = true;

    auto cov = estimate_fixed_covariance(model, params, fit, lopts, copts);

    std::cout << "theta_hat = " << fit.par[0] << "\n";
    std::cout << "theta_se  = " << params.params[0].std_error << "\n";
    std::cout << "cov(0,0)  = " << cov.covariance(0, 0) << "\n";

    // For this test model, profiled objective is:
    // 0.5 theta^2 + 0.5 n_random theta
    // Hessian = 1, covariance = 1, SE = 1.
    if (std::abs(params.params[0].std_error - 1.0) > 1e-3) {
        std::cerr << "FAIL: expected SE close to 1.0\n";
        return 1;
    }

    std::cout << "PASS\n";
    return 0;
}
