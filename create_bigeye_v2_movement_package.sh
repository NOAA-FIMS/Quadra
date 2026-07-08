#!/usr/bin/env bash
set -euo pipefail

BASE="examples/NMFS/pifsc_bigeye_tuna/v2"

mkdir -p "$BASE/architecture/steps/movement"
mkdir -p "$BASE/architecture/packages/movement"

cat > "$BASE/architecture/packages/movement/movement_context.hpp" <<'CPP'
#pragma once

#include "../../parameters/movement_parameters.hpp"
#include "../../state/population_state.hpp"

#include <vector>

namespace bigeye_v2 {

template <typename T>
struct MovementContext {
  const MovementParameters<T> *parameters = nullptr;
  std::vector<PopulationState<T>> *populations = nullptr;
};

} // namespace bigeye_v2
CPP

cat > "$BASE/architecture/steps/movement/identity_movement.hpp" <<'CPP'
#pragma once

#include "../../parameters/movement_parameters.hpp"
#include "../../state/population_state.hpp"

#include <vector>

namespace bigeye_v2 {

// Identity movement.
// Leaves population numbers unchanged.
struct IdentityMovement {
  template <typename T>
  void operator()(const MovementParameters<T> &,
                  std::vector<PopulationState<T>> &) const {
    // Intentionally empty.
  }
};

} // namespace bigeye_v2
CPP

cat > "$BASE/architecture/packages/movement/movement_package.hpp" <<'CPP'
#pragma once

#include "movement_context.hpp"
#include "../../steps/movement/identity_movement.hpp"
#include "../../../common/model_data.hpp"

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
CPP

mkdir -p "$BASE/16_movement_package_caa"

cat > "$BASE/16_movement_package_caa/bigeye_v2_16_movement_package_caa_check.cpp" <<'CPP'
#include "../architecture/packages/movement/movement_package.hpp"

#include <array>
#include <cstdlib>
#include <iostream>
#include <vector>

int main() {
  using namespace bigeye_v2;

  BigeyeModelData<double> data;
  data.n_years = 2;

  MovementParameters<double> parameters;

  std::vector<PopulationState<double>> populations(2);
  for (auto &pop : populations) {
    pop.numbers_at_age.assign(
        static_cast<std::size_t>(data.n_years),
        std::array<double, kAges>{});
  }

  populations[0].numbers_at_age[0][0] = 1000.0;
  populations[1].numbers_at_age[0][0] = 500.0;

  MovementContext<double> context{&parameters, &populations};

  MovementPackage{}(data, context);

  if (populations[0].numbers_at_age[0][0] != 1000.0 ||
      populations[1].numbers_at_age[0][0] != 500.0) {
    std::cerr << "FAIL: identity movement changed population numbers\n";
    return EXIT_FAILURE;
  }

  std::cout << "PASSED: Bigeye v2 CAA MovementPackage regression\n";
  return EXIT_SUCCESS;
}
CPP

cat > run_bigeye_v2_16_movement_package_caa_check.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples

c++ -std=c++17 -O3 \
  -I. \
  examples/NMFS/pifsc_bigeye_tuna/v2/16_movement_package_caa/bigeye_v2_16_movement_package_caa_check.cpp \
  -o build/examples/bigeye_v2_16_movement_package_caa_check

./build/examples/bigeye_v2_16_movement_package_caa_check
SH

chmod +x run_bigeye_v2_16_movement_package_caa_check.sh

echo "created CAA MovementPackage"
