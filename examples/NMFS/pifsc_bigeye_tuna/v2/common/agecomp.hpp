#pragma once

#include "bigeye_constants.hpp"
#include "derived.hpp"
#include "model_data.hpp"
#include "model_parameters.hpp"

namespace bigeye_v2 {

struct CatchAgeCompositionPrediction {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const BigeyeModelParameters<T> &, BigeyeDerived<T> &d) const {
    const auto ny = static_cast<std::size_t>(data.n_years);
    d.predicted_catch_age_proportion.assign(ny, std::array<T, kAges>{});

    for (std::size_t y = 0; y < ny; ++y) {
      T total = T(0.0);
      for (int a = 0; a < kAges; ++a) {
        total += d.catch_numbers_at_age[y][a];
      }

      for (int a = 0; a < kAges; ++a) {
        d.predicted_catch_age_proportion[y][a] =
            d.catch_numbers_at_age[y][a] / total;
      }
    }
  }
};

} // namespace bigeye_v2
