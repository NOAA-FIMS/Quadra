#include "test_common.hpp"
#include "../core/inference/covariance.hpp"
#include "../core/inference/report.hpp"

struct ReportV2Model {
    int n_random;

    template <typename T>
    T objective(const std::vector<T>& p) const {
        T theta = p[0];
        T nll = 0.5 * theta * theta;
        for (int i = 0; i < n_random; ++i) {
            T u = p[1 + i];
            nll += 0.5 * exp(theta) * u * u;
        }
        return nll;
    }

    template <typename T>
    T operator()(const std::vector<T>& p) const { return objective(p); }

    template <typename T, typename ReportLike>
    void report(const std::vector<T>& p, ReportLike& out) const {
        T theta = p[0];

        quadra::ReportMetadata dimless;
        dimless.units = "dimensionless";
        dimless.description = "Example scalar derived quantity";

        out.add("parameters/theta", theta, dimless);
        out.estimate("derived/exp_theta", exp(theta), dimless);
        out.estimate("derived/double_theta", T(2.0) * theta, dimless);
    }
};

int main() {
    using namespace quadra;

    quadra_tests::print_banner("Testing report v2 registry");

    int n_random = 10;
    ParameterVector params;
    params.add({"theta", 0.0, Transform::Identity, false});

    for (int i = 0; i < n_random; ++i) {
        params.add({"u_" + std::to_string(i), 0.0, Transform::Identity, true});
    }

    ReportV2Model model{n_random};

    LaplaceOptions lopts = quadra_tests::default_test_options();
    lopts.use_hutchinson_trace = false;

    auto fit = optimize_lbfgs(model, params, lopts);

    CovarianceOptions copts;
    copts.verbose = false;
    auto cov = estimate_fixed_covariance(model, params, fit, lopts, copts);

    auto report = evaluate_report_with_uncertainty(model, params, cov);
    report.print();
    report.to_csv("report_v2_test.csv");

    if (report.size() != 3) return 1;

    if (report.entries[1].metadata.group != "derived") {
        std::cerr << "FAIL: expected group derived\n";
        return 1;
    }

    double se_double_theta = report.entries[2].std_error;
    if (std::abs(se_double_theta - 2.0) > 1e-3) {
        std::cerr << "FAIL: expected SE(double_theta) close to 2\n";
        return 1;
    }

    std::cout << "PASS\n";
    return 0;
}
