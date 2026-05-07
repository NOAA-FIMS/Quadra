#include <cmath>
#include <iostream>
#include <vector>

#include "../core/optimizer/optimizer.hpp"
#include "../core/inference/covariance.hpp"
#include "../core/inference/report.hpp"

struct FisheriesReportV2Example {
    std::vector<double> log_index;
    std::vector<double> log_biomass;

    template <typename T>
    T objective(const std::vector<T>& p) const {
        T log_q = p[0];
        T log_sigma_obs = p[1];
        T log_sigma_re = p[2];

        T sigma_obs = exp(log_sigma_obs);
        T sigma_re = exp(log_sigma_re);

        T nll = 0.0;
        for (int t = 0; t < static_cast<int>(log_index.size()); ++t) {
            T u_t = p[3 + t];
            T pred = log_q + log_biomass[t] + u_t;
            T r = (log_index[t] - pred) / sigma_obs;

            nll += 0.5 * r * r + log_sigma_obs;

            T z = u_t / sigma_re;
            nll += 0.5 * z * z + log_sigma_re;
        }

        nll += 0.001 * log_q * log_q;
        nll += 0.001 * log_sigma_obs * log_sigma_obs;
        nll += 0.001 * log_sigma_re * log_sigma_re;
        return nll;
    }

    template <typename T>
    T operator()(const std::vector<T>& p) const { return objective(p); }

    template <typename T, typename ReportLike>
    void report(const std::vector<T>& p, ReportLike& out) const {
        quadra::ReportMetadata q_meta;
        q_meta.units = "index per biomass";
        q_meta.description = "Catchability coefficient";

        quadra::ReportMetadata index_meta;
        index_meta.units = "index units";
        index_meta.description = "Predicted index for first year";

        T q = exp(p[0]);

        out.add("parameters/log_q", p[0]);
        out.estimate("catchability/q", q, q_meta);
        out.estimate("predictions/index_year_0", q * exp(log_biomass[0] + p[3]), index_meta);
    }
};

int main() {
    using namespace quadra;

    std::vector<double> log_biomass = {
        std::log(1000), std::log(950), std::log(900), std::log(860)
    };
    std::vector<double> log_index = {
        std::log(10.1), std::log(9.4), std::log(9.2), std::log(8.5)
    };

    ParameterVector params;
    params.add({"log_q", std::log(0.01), Transform::Identity, false});
    params.add({"log_sigma_obs", std::log(0.2), Transform::Identity, false});
    params.add({"log_sigma_re", std::log(0.1), Transform::Identity, false});

    for (int t = 0; t < static_cast<int>(log_index.size()); ++t) {
        params.add({"u_year_" + std::to_string(t), 0.0, Transform::Identity, true});
    }

    FisheriesReportV2Example model{log_index, log_biomass};

    LaplaceOptions lopts;
    lopts.use_hutchinson_trace = false;
    lopts.hessian_drop_tol = 0.0;

    auto fit = optimize_lbfgs(model, params, lopts);

    CovarianceOptions copts;
    copts.verbose = false;
    auto cov = estimate_fixed_covariance(model, params, fit, lopts, copts);

    auto report = evaluate_report_with_uncertainty(model, params, cov);
    report.print();
    report.to_csv("fisheries_report.csv");

    return 0;
}
