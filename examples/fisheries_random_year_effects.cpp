
#include <cmath>
#include <iostream>
#include <vector>

#include "../core/optimizer.hpp"

DECLARE_ADGRAPH();

// Simple fisheries-style biomass index model with year random effects.
//
// log I_t ~ Normal(log(q) + log(B_t) + u_t, sigma_obs)
// u_t follows iid normal prior.
// This is not a full stock assessment model; it is a compact Laplace example.

struct FisheriesRandomYearEffects {
    std::vector<double> log_index;
    std::vector<double> log_biomass;

    template <typename T>
    T operator()(const std::vector<T>& p) const {
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

        // weak fixed-effect penalties
        nll += 0.001 * log_q * log_q;
        nll += 0.001 * log_sigma_obs * log_sigma_obs;
        nll += 0.001 * log_sigma_re * log_sigma_re;

        return nll;
    }
};

int main() {
    using namespace quadra;

    std::vector<double> log_biomass = {
        std::log(1000), std::log(950), std::log(900), std::log(860),
        std::log(830), std::log(790), std::log(760), std::log(720)
    };

    std::vector<double> log_index = {
        std::log(10.1), std::log(9.4), std::log(9.2), std::log(8.5),
        std::log(8.1), std::log(7.8), std::log(7.6), std::log(7.1)
    };

    ParameterVector params;
    params.add({"log_q", std::log(0.01), Transform::Identity, false});
    params.add({"log_sigma_obs", std::log(0.2), Transform::Identity, false});
    params.add({"log_sigma_re", std::log(0.1), Transform::Identity, false});

    for (int t = 0; t < static_cast<int>(log_index.size()); ++t) {
        params.add({"u_year_" + std::to_string(t), 0.0, Transform::Identity, true});
    }

    FisheriesRandomYearEffects model{log_index, log_biomass};

    LaplaceOptions opts;
    opts.use_hutchinson_trace = false; // small example: exact deterministic trace
    opts.hessian_drop_tol = 0.0;

    auto fit = optimize_lbfgs(model, params, opts);

    std::cout << "fit.value = " << fit.value << "\n";
    std::cout << "log_q_hat = " << fit.par[0] << "\n";

    return 0;
}
