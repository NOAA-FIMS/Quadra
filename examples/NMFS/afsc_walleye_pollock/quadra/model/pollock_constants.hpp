#pragma once

namespace pollock {

inline constexpr int n_ages = 7;

inline constexpr double weight_at_age[n_ages] = {
    0.20, 0.45, 0.75, 1.10, 1.45, 1.75, 2.00};

inline constexpr double maturity_at_age[n_ages] = {
    0.00, 0.10, 0.45, 0.80, 0.95, 1.00, 1.00};

inline constexpr double natural_mortality = 0.25;

}  // namespace pollock
