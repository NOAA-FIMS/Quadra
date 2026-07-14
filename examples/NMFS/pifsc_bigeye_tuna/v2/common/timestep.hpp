#pragma once

#include "bigeye_constants.hpp"
#include "fleet.hpp"
#include "model_data.hpp"

#include <array>
#include <cmath>
#include <vector>

namespace bigeye_v2 {

template <typename T> struct TimeStepPopulationState {
  std::vector<std::array<T, kAges>> numbers_at_age{};
};

struct AnnualTimeStep {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data, const BigeyeDerived<T> &life,
                  const BigeyeFleetDerived<T> &fleet, const T &recruitment,
                  TimeStepPopulationState<T> &pop) const {
    const auto ny = static_cast<std::size_t>(data.n_years);

    if (pop.numbers_at_age.size() < ny) {
      pop.numbers_at_age.assign(ny, std::array<T, kAges>{});
    }

    for (std::size_t y = 0; y + 1 < ny; ++y) {
      pop.numbers_at_age[y + 1][0] = recruitment;

      for (int a = 1; a < kAges - 1; ++a) {
        pop.numbers_at_age[y + 1][a] =
            pop.numbers_at_age[y][a - 1] * std::exp(-fleet.z_at_age[a - 1]);
      }

      const int plus = kAges - 1;

      pop.numbers_at_age[y + 1][plus] =
          pop.numbers_at_age[y][plus - 1] *
              std::exp(-fleet.z_at_age[plus - 1]) +
          pop.numbers_at_age[y][plus] * std::exp(-fleet.z_at_age[plus]);
    }
  }
};

struct AnnualCatchAtStartOfYear {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data, const BigeyeDerived<T> &life,
                  const BigeyeFleetDerived<T> &fleet_state,
                  const TimeStepPopulationState<T> &pop,
                  BigeyeFleetDerived<T> &fleet) const {
    const auto ny = static_cast<std::size_t>(data.n_years);

    fleet.catch_numbers_at_age.assign(ny, std::array<T, kAges>{});
    fleet.catch_biomass_at_age.assign(ny, std::array<T, kAges>{});
    fleet.total_catch_biomass_by_year.assign(ny, T(0.0));

    for (std::size_t y = 0; y < ny; ++y) {
      for (int a = 0; a < kAges; ++a) {
        const T z = fleet_state.z_at_age[a];

        const T cn = pop.numbers_at_age[y][a] * (fleet_state.f_at_age[a] / z) *
                     (T(1.0) - std::exp(-z));

        fleet.catch_numbers_at_age[y][a] = cn;
        fleet.catch_biomass_at_age[y][a] = cn * life.weight_at_age[a];
        fleet.total_catch_biomass_by_year[y] +=
            fleet.catch_biomass_at_age[y][a];
      }
    }
  }
};

} // namespace bigeye_v2
