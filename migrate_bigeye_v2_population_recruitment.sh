#!/usr/bin/env bash
set -euo pipefail

BASE="examples/NMFS/pifsc_bigeye_tuna/v2"

cat > "$BASE/architecture/steps/population/recruitment.hpp" <<'CPP'
#pragma once

#include "../../parameters/population_parameters.hpp"
#include "../../state/population_state.hpp"
#include "../../../common/model_data.hpp"

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
    population.recruits_by_year.assign(
        static_cast<std::size_t>(data.n_years), p.r0);
  }
};

} // namespace bigeye_v2
CPP

mkdir -p "$BASE/04_population_recruitment_caa"

cat > "$BASE/04_population_recruitment_caa/bigeye_v2_04_population_recruitment_caa_check.cpp" <<'CPP'
#include "../architecture/parameters/population_parameters.hpp"
#include "../architecture/state/population_state.hpp"
#include "../architecture/steps/population/recruitment.hpp"
#include "../common/model_data.hpp"

#include <cstdlib>
#include <iostream>

int main() {
  using namespace bigeye_v2;

  BigeyeModelData<double> data;
  data.n_years = 4;

  PopulationParameters<double> p;
  p.r0 = 1234.5;

  PopulationState<double> population;

  FixedRecruitment{}(data, p, population);

  if (population.recruits_by_year.size() != 4) {
    std::cerr << "FAIL: recruits size\n";
    return EXIT_FAILURE;
  }

  for (double r : population.recruits_by_year) {
    if (r != 1234.5) {
      std::cerr << "FAIL: recruit value\n";
      return EXIT_FAILURE;
    }
  }

  std::cout << "PASSED: Bigeye v2 CAA population recruitment regression\n";
  return EXIT_SUCCESS;
}
CPP

cat > run_bigeye_v2_04_population_recruitment_caa_check.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples

c++ -std=c++17 -O3 \
  -I. \
  examples/NMFS/pifsc_bigeye_tuna/v2/04_population_recruitment_caa/bigeye_v2_04_population_recruitment_caa_check.cpp \
  -o build/examples/bigeye_v2_04_population_recruitment_caa_check

./build/examples/bigeye_v2_04_population_recruitment_caa_check
SH

chmod +x run_bigeye_v2_04_population_recruitment_caa_check.sh

echo "migrated CAA population recruitment step"
