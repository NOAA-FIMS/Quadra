#pragma once

#include "../../../common/model_data.hpp"
#include "../../parameters/population_parameters.hpp"
#include "../../state/population_state.hpp"

//------------------------------------------------------------
// Recruitment
//
// Purpose
// -------
// Fills annual recruits for one population.
//
// Consumes
// --------
// BigeyeModelData
// PopulationParameters
//
// Produces
// --------
// PopulationState::recruits_by_year
//
// Notes
// -----
// Stateless.
// Owns no memory beyond resizing population-owned state.
//------------------------------------------------------------

namespace bigeye_v2 {

struct FixedRecruitment {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const PopulationParameters<T> &p,
                  PopulationState<T> &population) const {
    population.recruits_by_year.assign(static_cast<std::size_t>(data.n_years),
                                       p.r0);
  }
};

} // namespace bigeye_v2
