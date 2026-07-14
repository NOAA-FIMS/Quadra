#pragma once

#include "../../../common/model_data.hpp"
#include "../../parameters/population_parameters.hpp"
#include "../../state/fleet_state.hpp"
#include "../../state/life_history_state.hpp"
#include "../../state/population_state.hpp"
#include "../../steps/population/aging.hpp"
#include "../../steps/population/plus_group.hpp"
#include "../../steps/population/recruitment.hpp"
#include "../../steps/population/spawning_biomass.hpp"
#include "../../steps/population/survival.hpp"
#include "population_context.hpp"

#include <array>

namespace bigeye_v2 {

//------------------------------------------------------------
// PopulationPackage
//
// Purpose
// -------
// Composes the population steps for one population.
//
// Sequence
// --------
// Recruit
// Survive
// Age
// Accumulate plus group
// Compute spawning biomass
//
// Consumes
// --------
// BigeyeModelData
// PopulationParameters
// LifeHistoryState
// FleetState::z_at_age
//
// Produces
// --------
// PopulationState::recruits_by_year
// PopulationState::survivors_at_age
// PopulationState::numbers_at_age
// PopulationState::spawning_biomass_by_year
//
// Notes
// -----
// State owns memory.
// Steps own algorithms.
// Packages orchestrate related steps.
//------------------------------------------------------------
struct PopulationPackage {

  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const PopulationContext<T> &context) const {
    (*this)(data, *context.parameters, *context.life, *context.fleet,
            *context.population);
  }

  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const PopulationParameters<T> &parameters,
                  const LifeHistoryState<T> &life, const FleetState<T> &fleet,
                  PopulationState<T> &population) const {
    FixedRecruitment{}(data, parameters, population);

    if (population.numbers_at_age.size() <
        static_cast<std::size_t>(data.n_years)) {
      population.numbers_at_age.assign(static_cast<std::size_t>(data.n_years),
                                       std::array<T, kAges>{});
    }

    for (std::size_t y = 0; y < population.recruits_by_year.size(); ++y) {
      population.numbers_at_age[y][0] = population.recruits_by_year[y];
    }

    Survival{}(fleet, population);
    Aging{}(population);

    for (std::size_t y = 0; y < population.recruits_by_year.size(); ++y) {
      population.numbers_at_age[y][0] = population.recruits_by_year[y];
    }

    PlusGroup{}(population);
    SpawningBiomass{}(life, population);
  }
};

} // namespace bigeye_v2
