#ifndef QUADRA_BIG_CATCH_AT_AGE_DERIVED_HPP
#define QUADRA_BIG_CATCH_AT_AGE_DERIVED_HPP

#include "catch_at_age_shared.hpp"

#include <cmath>
#include <limits>
#include <vector>

namespace example {

struct CatchAtAgeDerivedQuantities {
    double terminal_ssb_proxy_m = std::numeric_limits<double>::quiet_NaN();
    double terminal_depletion_m = std::numeric_limits<double>::quiet_NaN();
    double mean_f_m = std::numeric_limits<double>::quiet_NaN();
};

inline CatchAtAgeDerivedQuantities evaluate_catch_at_age_derived_quantities(
    const CatchAtAgeLaplaceModel& model,
    const std::vector<double>& theta,
    const std::vector<double>& u)
{
    CatchAtAgeDerivedQuantities out;

    if (theta.size() < 9 ||
        u.size() < static_cast<std::size_t>(model.data.n_years)) {
        return out;
    }

    const int n_years = model.data.n_years;
    const int n_ages = model.data.n_ages;

    const double R0 = std::exp(theta[0]);
    const double M = std::exp(theta[1]);
    const double Fbar = std::exp(theta[3]);
    const double sel50_raw = theta[4];
    const double sel_slope = std::exp(theta[5]);

    std::vector<double> selectivity(static_cast<std::size_t>(n_ages), 0.0);

    for (int a = 0; a < n_ages; ++a) {
        const double age = static_cast<double>(a + 1);
        const double sel50 = 3.0 + sel50_raw;

        selectivity[static_cast<std::size_t>(a)] =
            1.0 / (1.0 + std::exp(-(age - sel50) / sel_slope));
    }

    std::vector<double> numbers(
        static_cast<std::size_t>(n_ages),
        R0 / static_cast<double>(n_ages));

    double ssb0 = std::numeric_limits<double>::quiet_NaN();
    double terminal_ssb = std::numeric_limits<double>::quiet_NaN();

    double f_sum = 0.0;
    int f_count = 0;

    for (int y = 0; y < n_years; ++y) {
        const double rec = R0 * std::exp(u[static_cast<std::size_t>(y)]);
        numbers[0] = rec;

        double ssb = 0.0;

        for (int a = 0; a < n_ages; ++a) {
            const double age = static_cast<double>(a + 1);
            const double weight = (age * age * age) / 100.0;
            const double maturity =
                1.0 / (1.0 + std::exp(-(age - 4.0)));

            const double Na = numbers[static_cast<std::size_t>(a)];
            const double Sa = selectivity[static_cast<std::size_t>(a)];
            const double F_at_age = Fbar * Sa;

            ssb += Na * weight * maturity;
            f_sum += F_at_age;
            f_count++;
        }

        if (y == 0) {
            ssb0 = ssb;
        }

        terminal_ssb = ssb;

        for (int a = n_ages - 1; a >= 1; --a) {
            const double Sa_prev = selectivity[static_cast<std::size_t>(a - 1)];
            const double Z_prev = M + Fbar * Sa_prev;

            numbers[static_cast<std::size_t>(a)] =
                numbers[static_cast<std::size_t>(a - 1)] * std::exp(-Z_prev);
        }
    }

    out.terminal_ssb_proxy_m = terminal_ssb;

    if (std::isfinite(ssb0) && ssb0 > 0.0) {
        out.terminal_depletion_m = terminal_ssb / ssb0;
    }

    if (f_count > 0) {
        out.mean_f_m = f_sum / static_cast<double>(f_count);
    }

    return out;
}

} // namespace example

#endif
