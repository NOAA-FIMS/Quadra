#pragma once

#include "../../parameters/fleet_parameters.hpp"
#include "../../state/fleet_state.hpp"
#include "../../state/likelihood_state.hpp"

namespace bigeye_v2 {

template <typename T>
struct LikelihoodContext {
  const FleetParameters<T> *parameters = nullptr;
  const FleetState<T> *fleet = nullptr;
  LikelihoodState<T> *likelihood = nullptr;
};

} // namespace bigeye_v2
