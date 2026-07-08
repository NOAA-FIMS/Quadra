#pragma once

#include "derived.hpp"
#include "model_data.hpp"
#include "model_parameters.hpp"

#include <vector>

namespace bigeye_v2 {

struct FixedRecruitment {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const BigeyeModelParameters<T> &p,
                  BigeyeDerived<T> &d) const {
    d.recruits_by_year.assign(static_cast<std::size_t>(data.n_years), p.r0);
  }
};

} // namespace bigeye_v2
