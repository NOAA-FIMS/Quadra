
#include <cmath>
#include <iostream>
#include <vector>

#include "../core/optimizer.hpp"

// Fisheries-style selectivity example.
//
// Observed log selectivity-at-age is modeled with random-walk deviations.
// This is a compact example for sparse random-effect Hessians over age.

struct SelectivityRandomWalk {
    std::vector<double> obs_log_sel;

    template <typename T>
    T operator()(const std::vector<T>& p) const {
        T a50 = p[0];
        T log_slope = p[1];
        T log_sigma_obs = p[2];
        T log_sigma_rw = p[3];

        T slope = exp(log_slope);
        T sigma_obs = exp(log_sigma_obs);
        T sigma_rw = exp(log_sigma_rw);

        T nll = 0.0;

        for (int a = 0; a < static_cast<int>(obs_log_sel.size()); ++a) {
            T age = T(a + 1);
            T dev = p[4 + a];

            T sel = 1.0 / (1.0 + exp(-slope * (age - a50)));
            T log_sel = log(sel) + dev;

            T r = (obs_log_sel[a] - log_sel) / sigma_obs;
            nll += 0.5 * r * r + log_sigma_obs;

            if (a == 0) {
                T z = dev / sigma_rw;
                nll += 0.5 * z * z + log_sigma_rw;
            } else {
                T dev_prev = p[4 + a - 1];
                T z = (dev - dev_prev) / sigma_rw;
                nll += 0.5 * z * z + log_sigma_rw;
            }
        }

        nll += 0.001 * a50 * a50;
        nll += 0.001 * log_slope * log_slope;
        nll += 0.001 * log_sigma_obs * log_sigma_obs;
        nll += 0.001 * log_sigma_rw * log_sigma_rw;

        return nll;
    }
};

int main() {
    using namespace quadra;

    std::vector<double> obs_log_sel = {
        std::log(0.05), std::log(0.12), std::log(0.28),
        std::log(0.55), std::log(0.78), std::log(0.90),
        std::log(0.96), std::log(0.98)
    };

    ParameterVector params;
    params.add({"a50", 4.0, Transform::Identity, false});
    params.add({"log_slope", std::log(1.0), Transform::Identity, false});
    params.add({"log_sigma_obs", std::log(0.2), Transform::Identity, false});
    params.add({"log_sigma_rw", std::log(0.1), Transform::Identity, false});

    for (int a = 0; a < static_cast<int>(obs_log_sel.size()); ++a) {
        params.add({"sel_dev_age_" + std::to_string(a + 1), 0.0, Transform::Identity, true});
    }

    SelectivityRandomWalk model{obs_log_sel};

    LaplaceOptions opts;
    opts.use_hutchinson_trace = false;
    opts.hessian_drop_tol = 0.0;

    auto fit = optimize_lbfgs(model, params, opts);

    std::cout << "fit.value = " << fit.value << "\n";
    std::cout << "a50_hat = " << fit.par[0] << "\n";

    return 0;
}
