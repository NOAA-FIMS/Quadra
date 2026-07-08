#pragma once

#include "../../state/life_history_state.hpp"
#include "../../state/population_state.hpp"

namespace bigeye_v2 {

// Computes spawning biomass by year.
struct SpawningBiomass {
  template <typename T>
  void operator()(const LifeHistoryState<T> &life,
                  PopulationState<T> &population) const {
    const auto ny = population.numbers_at_age.size();

    population.spawning_biomass_by_year.assign(ny, T(0.0));

    for (std::size_t y = 0; y < ny; ++y) {
      for (int a = 0; a < kAges; ++a) {
        population.spawning_biomass_by_year[y] +=
            population.numbers_at_age[y][a] *
            life.weight_at_age[a] *
            life.maturity_at_age[a];
      }
    }
  }
};

} // namespace bigeye_v2
