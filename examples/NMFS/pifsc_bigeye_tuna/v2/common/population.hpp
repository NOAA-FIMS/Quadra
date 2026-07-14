#pragma once

#include "bigeye_constants.hpp"
#include "derived.hpp"
#include "model_data.hpp"
#include "model_parameters.hpp"

#include <cmath>
#include <vector>

namespace bigeye_v2 {

struct UnfishedPopulation {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const BigeyeModelParameters<T> & /*p*/,
                  BigeyeDerived<T> &d) const {
    const auto ny = static_cast<std::size_t>(data.n_years);

    d.numbers_at_age.assign(ny, std::array<T, kAges>{});
    d.spawning_biomass_by_year.assign(ny, T(0.0));

    for (std::size_t y = 0; y < ny; ++y) {
      d.numbers_at_age[y][0] = d.recruits_by_year[y];

      for (int a = 1; a < kAges - 1; ++a) {
        d.numbers_at_age[y][a] =
            d.numbers_at_age[y][a - 1] * std::exp(-d.m_at_age[a - 1]);
      }

      const int plus = kAges - 1;
      d.numbers_at_age[y][plus] = d.numbers_at_age[y][plus - 1] *
                                  std::exp(-d.m_at_age[plus - 1]) /
                                  (T(1.0) - std::exp(-d.m_at_age[plus]));

      for (int a = 0; a < kAges; ++a) {
        d.spawning_biomass_by_year[y] +=
            d.numbers_at_age[y][a] * d.weight_at_age[a] * d.maturity_at_age[a];
      }
    }
  }
};

} // namespace bigeye_v2
