#pragma once

#include "../../../common/model_data.hpp"
#include "../../steps/movement/identity_movement.hpp"
#include "movement_context.hpp"

namespace bigeye_v2 {

//------------------------------------------------------------
// MovementPackage
//
// Purpose
// -------
// Composes movement steps across populations.
//
// Sequence
// --------
// Move
//
// Notes
// -----
// State owns memory.
// Steps own algorithms.
// Packages orchestrate related steps.
//------------------------------------------------------------
struct MovementPackage {
  template <typename T>
  void operator()(const BigeyeModelData<T> &,
                  const MovementContext<T> &context) const {
    IdentityMovement{}(*context.parameters, *context.populations);
  }
};

} // namespace bigeye_v2
