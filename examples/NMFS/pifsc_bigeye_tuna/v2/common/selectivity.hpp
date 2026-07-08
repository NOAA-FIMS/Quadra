#pragma once

#include "bigeye_constants.hpp"
#include "derived.hpp"
#include "model_data.hpp"
#include "model_parameters.hpp"

#include <cmath>

namespace bigeye_v2 {

struct LogisticSelectivity {
  template <typename T>
  void operator()(const BigeyeModelData<T> & /*data*/,
                  const BigeyeModelParameters<T> &p,
                  BigeyeDerived<T> &d) const {
    for (int a = 0; a < kAges; ++a) {
      const T age = T(a + 1);
      d.selectivity_at_age[a] =
          T(1.0) / (T(1.0) + std::exp(-p.sel_slope * (age - p.sel_a50)));
    }
  }
};

} // namespace bigeye_v2
