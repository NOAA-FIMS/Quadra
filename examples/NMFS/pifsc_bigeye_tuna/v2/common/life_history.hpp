#pragma once

#include "bigeye_constants.hpp"
#include "derived.hpp"
#include "model_data.hpp"
#include "model_parameters.hpp"

#include <cmath>

namespace bigeye_v2 {

struct BigeyeLifeHistory {
  template <typename T>
  void operator()(const BigeyeModelData<T> & /*data*/,
                  const BigeyeModelParameters<T> &p,
                  BigeyeDerived<T> &d) const {
    const T adult_m = T(0.45);
    const T m_young = adult_m * std::exp(p.log_m_young_offset);
    const T m_old = adult_m * std::exp(p.log_m_old_offset);

    static constexpr double weight_values[kAges] = {
        1.0, 3.0, 5.0, 7.0, 9.0, 11.0, 13.0, 15.0, 17.0, 19.0};

    static constexpr double maturity_values[kAges] = {
        0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};

    for (int a = 0; a < kAges; ++a) {
      if (a <= 2) {
        d.m_at_age[a] = m_young;
      } else if (a >= 7) {
        d.m_at_age[a] = m_old;
      } else {
        d.m_at_age[a] = adult_m;
      }

      d.weight_at_age[a] = T(weight_values[a]);
      d.maturity_at_age[a] = T(maturity_values[a]);
    }
  }
};

} // namespace bigeye_v2
