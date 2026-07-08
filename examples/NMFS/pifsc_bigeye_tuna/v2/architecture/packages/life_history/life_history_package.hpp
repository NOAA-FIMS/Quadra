#pragma once

#include "../../parameters/life_history_parameters.hpp"
#include "../../state/life_history_state.hpp"
#include "../../steps/life_history/life_history.hpp"
#include "../../../common/model_data.hpp"

namespace bigeye_v2 {

struct LifeHistoryPackage {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const LifeHistoryParameters<T> &parameters,
                  LifeHistoryState<T> &life) const {
    BigeyeLifeHistory{}(data, parameters, life);
  }
};

} // namespace bigeye_v2
