#include "../examples/big/catch_at_age_derived.hpp"

#include <cmath>
#include <iostream>
#include <vector>

int main()
{
    example::CatchAtAgeLaplaceModel model;

    std::vector<double> theta{
        std::log(900.0),
        std::log(0.25),
        std::log(0.15),
        std::log(0.18),
        0.0,
        std::log(1.25),
        std::log(0.20),
        std::log(0.18),
        std::log(0.35)
    };

    std::vector<double> u(static_cast<std::size_t>(model.data.n_years), 0.0);

    auto derived =
        example::evaluate_catch_at_age_derived_quantities(
            model,
            theta,
            u);

    if (!std::isfinite(derived.terminal_ssb_proxy_m) ||
        !std::isfinite(derived.terminal_depletion_m) ||
        !std::isfinite(derived.mean_f_m)) {
        std::cerr << "FAIL: non-finite derived quantity\n";
        return 1;
    }

    if (!(derived.terminal_ssb_proxy_m > 0.0)) {
        std::cerr << "FAIL: terminal SSB proxy is not positive\n";
        return 1;
    }

    if (!(derived.terminal_depletion_m > 0.0)) {
        std::cerr << "FAIL: terminal depletion is not positive\n";
        return 1;
    }

    if (!(derived.mean_f_m > 0.0)) {
        std::cerr << "FAIL: mean F is not positive\n";
        return 1;
    }

    std::cout << "PASS: big catch-at-age derived quantities\n";
    std::cout << "  terminal_ssb_proxy: "
              << derived.terminal_ssb_proxy_m << "\n";
    std::cout << "  terminal_depletion: "
              << derived.terminal_depletion_m << "\n";
    std::cout << "  mean_f: "
              << derived.mean_f_m << "\n";

    return 0;
}
