
#include "test_common.hpp"

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
    quadra_tests::print_banner("Testing exact Hdot validation path");

#ifndef QUADRA_VALIDATE_HDOT
    std::cout << "NOTE: Compile with -DQUADRA_VALIDATE_HDOT to print exact-vs-FD Hdot checks.\n";
#endif

    int n_random = 10;

    ParameterVector params;
    params.add({"theta", 0.0, Transform::Identity, false});

    for (int i = 0; i < n_random; ++i) {
        params.add({"u_" + std::to_string(i), 0.0, Transform::Identity, true});
    }

    CurvatureDependsOnTheta model{n_random};

    auto opts = quadra_tests::default_test_options();
    opts.validate_hdot = true;

    auto fit = optimize_lbfgs(model, params, opts);

    std::cout << "fit.value = " << fit.value << "\n";
    std::cout << "PASS\n";
    return 0;
}
