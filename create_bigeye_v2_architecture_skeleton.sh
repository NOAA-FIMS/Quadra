#!/usr/bin/env bash
set -euo pipefail

BASE="examples/NMFS/pifsc_bigeye_tuna/v2"
mkdir -p "$BASE/architecture/state" \
         "$BASE/architecture/parameters" \
         "$BASE/architecture/operators"

cat > "$BASE/architecture/state/population_state.hpp" <<'CPP'
#pragma once
#include "../../common/bigeye_constants.hpp"
#include <array>
#include <vector>

namespace bigeye_v2 {
template <typename T>
struct PopulationState {
  std::vector<T> recruits_by_year{};
  std::vector<std::array<T, kAges>> numbers_at_age{};
  std::vector<T> spawning_biomass_by_year{};
};
}
CPP

cat > "$BASE/architecture/state/fleet_state.hpp" <<'CPP'
#pragma once
#include "../../common/bigeye_constants.hpp"
#include <array>
#include <vector>

namespace bigeye_v2 {
template <typename T>
struct FleetState {
  std::array<T, kAges> selectivity_at_age{};
  std::array<T, kAges> f_at_age{};
  std::array<T, kAges> z_at_age{};
  std::vector<std::array<T, kAges>> catch_numbers_at_age{};
  std::vector<std::array<T, kAges>> catch_biomass_at_age{};
  std::vector<T> total_catch_biomass_by_year{};
};
}
CPP

cat > "$BASE/architecture/parameters/fleet_parameters.hpp" <<'CPP'
#pragma once
namespace bigeye_v2 {
template <typename T>
struct FleetParameters {
  T sel_a50 = T(5.0);
  T sel_slope = T(1.0);
  T fbar = T(0.2);
};
}
CPP

cat > "$BASE/architecture/operators/fleet_operators.hpp" <<'CPP'
#pragma once
#include "../../common/bigeye_constants.hpp"
#include "../../common/derived.hpp"
#include "../../common/model_data.hpp"
#include "../parameters/fleet_parameters.hpp"
#include "../state/fleet_state.hpp"
#include "../state/population_state.hpp"
#include <cmath>

namespace bigeye_v2 {

struct FleetLogisticSelectivity {
  template <typename T>
  void operator()(const BigeyeModelData<T>&,
                  const FleetParameters<T>& p,
                  FleetState<T>& fleet) const {
    for (int a = 0; a < kAges; ++a) {
      const T age = T(a + 1);
      fleet.selectivity_at_age[a] =
          T(1.0) / (T(1.0) + std::exp(-p.sel_slope * (age - p.sel_a50)));
    }
  }
};

struct FleetFishingMortality {
  template <typename T>
  void operator()(const BigeyeModelData<T>&,
                  const FleetParameters<T>& p,
                  const BigeyeDerived<T>& life,
                  FleetState<T>& fleet) const {
    for (int a = 0; a < kAges; ++a) {
      fleet.f_at_age[a] = p.fbar * fleet.selectivity_at_age[a];
      fleet.z_at_age[a] = life.m_at_age[a] + fleet.f_at_age[a];
    }
  }
};

}
CPP

echo "created Bigeye v2 architecture skeleton"
