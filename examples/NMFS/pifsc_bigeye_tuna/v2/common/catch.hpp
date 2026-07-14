#pragma once

#include "bigeye_constants.hpp"
#include "derived.hpp"
#include "model_data.hpp"
#include "model_parameters.hpp"

#include <array>
#include <cmath>
#include <vector>

namespace bigeye_v2 {

struct BaranovCatch {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const BigeyeModelParameters<T> & /*p*/,
                  BigeyeDerived<T> &d) const {
    const auto ny = static_cast<std::size_t>(data.n_years);

    d.catch_numbers_at_age.assign(ny, std::array<T, kAges>{});
    d.catch_biomass_at_age.assign(ny, std::array<T, kAges>{});
    d.total_catch_biomass_by_year.assign(ny, T(0.0));

    for (std::size_t y = 0; y < ny; ++y) {
      for (int a = 0; a < kAges; ++a) {
        const T z = d.z_at_age[a];

        const T catch_numbers = d.numbers_at_age[y][a] * (d.f_at_age[a] / z) *
                                (T(1.0) - std::exp(-z));

        d.catch_numbers_at_age[y][a] = catch_numbers;
        d.catch_biomass_at_age[y][a] = catch_numbers * d.weight_at_age[a];
        d.total_catch_biomass_by_year[y] += d.catch_biomass_at_age[y][a];
      }
    }
  }
};

} // namespace bigeye_v2
