#!/usr/bin/env bash
set -euo pipefail

BASE="examples/NMFS/pifsc_bigeye_tuna/v2"

mkdir -p "$BASE/architecture/steps/likelihood"
mkdir -p "$BASE/architecture/packages/likelihood"

cat > "$BASE/architecture/steps/likelihood/lognormal_catch.hpp" <<'CPP'
#pragma once

#include "../../parameters/fleet_parameters.hpp"
#include "../../state/fleet_state.hpp"
#include "../../state/likelihood_state.hpp"
#include "../../../common/model_data.hpp"

#include <cmath>

namespace bigeye_v2 {

struct LognormalCatchLikelihood {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const FleetParameters<T> &parameters,
                  const FleetState<T> &fleet,
                  LikelihoodState<T> &likelihood) const {
    likelihood.catch_nll = T(0.0);

    for (std::size_t y = 0; y < data.observed_catch_biomass_by_year.size(); ++y) {
      const T obs = data.observed_catch_biomass_by_year[y];
      const T pred = fleet.total_catch_biomass_by_year[y];
      const T r = (std::log(obs) - std::log(pred)) / parameters.catch_sigma;

      likelihood.catch_nll +=
          T(0.5) * r * r +
          std::log(parameters.catch_sigma) +
          std::log(obs);
    }

    likelihood.total_nll += likelihood.catch_nll;
  }
};

} // namespace bigeye_v2
CPP

cat > "$BASE/architecture/steps/likelihood/lognormal_index.hpp" <<'CPP'
#pragma once

#include "../../parameters/fleet_parameters.hpp"
#include "../../state/fleet_state.hpp"
#include "../../state/likelihood_state.hpp"
#include "../../../common/model_data.hpp"

#include <cmath>

namespace bigeye_v2 {

struct LognormalIndexLikelihood {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const FleetParameters<T> &parameters,
                  const FleetState<T> &fleet,
                  LikelihoodState<T> &likelihood) const {
    likelihood.index_nll = T(0.0);

    for (std::size_t y = 0; y < data.observed_index_by_year.size(); ++y) {
      const T obs = data.observed_index_by_year[y];
      const T pred = fleet.predicted_index_by_year[y];
      const T r = (std::log(obs) - std::log(pred)) / parameters.index_sigma;

      likelihood.index_nll +=
          T(0.5) * r * r +
          std::log(parameters.index_sigma) +
          std::log(obs);
    }

    likelihood.total_nll += likelihood.index_nll;
  }
};

} // namespace bigeye_v2
CPP

cat > "$BASE/architecture/steps/likelihood/multinomial_agecomp.hpp" <<'CPP'
#pragma once

#include "../../state/fleet_state.hpp"
#include "../../state/likelihood_state.hpp"
#include "../../../common/bigeye_constants.hpp"
#include "../../../common/model_data.hpp"

#include <cmath>

namespace bigeye_v2 {

struct MultinomialAgeCompLikelihood {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const FleetState<T> &fleet,
                  LikelihoodState<T> &likelihood) const {
    likelihood.agecomp_nll = T(0.0);

    const T eps = T(1.0e-12);

    for (std::size_t y = 0; y < data.observed_catch_age_proportion.size(); ++y) {
      const T n_eff = data.catch_agecomp_sample_size[y];

      for (int a = 0; a < kAges; ++a) {
        const T obs = data.observed_catch_age_proportion[y][a];
        const T pred = fleet.predicted_catch_age_proportion[y][a] + eps;

        likelihood.agecomp_nll -= n_eff * obs * std::log(pred);
      }
    }

    likelihood.total_nll += likelihood.agecomp_nll;
  }
};

} // namespace bigeye_v2
CPP

cat > "$BASE/architecture/packages/likelihood/likelihood_package.hpp" <<'CPP'
#pragma once

#include "../../parameters/fleet_parameters.hpp"
#include "../../state/fleet_state.hpp"
#include "../../state/likelihood_state.hpp"
#include "../../steps/likelihood/lognormal_catch.hpp"
#include "../../steps/likelihood/lognormal_index.hpp"
#include "../../steps/likelihood/multinomial_agecomp.hpp"
#include "../../../common/model_data.hpp"

namespace bigeye_v2 {

struct LikelihoodPackage {
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const FleetParameters<T> &parameters,
                  const FleetState<T> &fleet,
                  LikelihoodState<T> &likelihood) const {
    likelihood.catch_nll = T(0.0);
    likelihood.index_nll = T(0.0);
    likelihood.agecomp_nll = T(0.0);
    likelihood.total_nll = T(0.0);

    LognormalCatchLikelihood{}(data, parameters, fleet, likelihood);
    LognormalIndexLikelihood{}(data, parameters, fleet, likelihood);
    MultinomialAgeCompLikelihood{}(data, fleet, likelihood);
  }
};

} // namespace bigeye_v2
CPP

mkdir -p "$BASE/14_likelihood_package_caa"

cat > "$BASE/14_likelihood_package_caa/bigeye_v2_14_likelihood_package_caa_check.cpp" <<'CPP'
#include "../architecture/packages/likelihood/likelihood_package.hpp"

#include <array>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>

namespace {
bool nearly_equal(double a, double b, double tol = 1.0e-6) {
  return std::abs(a - b) <= tol;
}
} // namespace

int main() {
  using namespace bigeye_v2;

  BigeyeModelData<double> data;
  data.n_years = 3;

  data.observed_catch_biomass_by_year = {100.0, 120.0, 90.0};
  data.observed_index_by_year = {10.0, 12.0, 9.0};
  data.catch_agecomp_sample_size = {100.0, 100.0, 100.0};

  data.observed_catch_age_proportion = {
      std::array<double, kAges>{0.02, 0.04, 0.07, 0.10, 0.14,
                                0.16, 0.17, 0.14, 0.10, 0.06},
      std::array<double, kAges>{0.02, 0.04, 0.07, 0.10, 0.14,
                                0.16, 0.17, 0.14, 0.10, 0.06},
      std::array<double, kAges>{0.02, 0.04, 0.07, 0.10, 0.14,
                                0.16, 0.17, 0.14, 0.10, 0.06}};

  FleetParameters<double> parameters;
  parameters.catch_sigma = 0.1;
  parameters.index_sigma = 0.2;

  FleetState<double> fleet;
  fleet.total_catch_biomass_by_year = {100.0, 120.0, 90.0};
  fleet.predicted_index_by_year = {10.0, 12.0, 9.0};

  fleet.predicted_catch_age_proportion.assign(
      static_cast<std::size_t>(data.n_years),
      std::array<double, kAges>{});

  for (int y = 0; y < data.n_years; ++y) {
    fleet.predicted_catch_age_proportion[y] =
        std::array<double, kAges>{0.02, 0.04, 0.07, 0.10, 0.14,
                                  0.16, 0.17, 0.14, 0.10, 0.06};
  }

  LikelihoodState<double> likelihood;

  LikelihoodPackage{}(data, parameters, fleet, likelihood);

  constexpr double expected_catch_nll = 7.783640596221253;
  constexpr double expected_index_nll = 5.886104031450156;
  constexpr double expected_agecomp_nll = 612.2363001098417;
  constexpr double expected_total_nll =
      expected_catch_nll + expected_index_nll + expected_agecomp_nll;

  if (!nearly_equal(likelihood.catch_nll, expected_catch_nll)) {
    std::cerr << std::setprecision(17)
              << "FAIL: catch_nll got " << likelihood.catch_nll
              << " expected " << expected_catch_nll
              << " diff " << (likelihood.catch_nll - expected_catch_nll)
              << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(likelihood.index_nll, expected_index_nll)) {
    std::cerr << std::setprecision(17)
              << "FAIL: index_nll got " << likelihood.index_nll
              << " expected " << expected_index_nll
              << " diff " << (likelihood.index_nll - expected_index_nll)
              << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(likelihood.agecomp_nll, expected_agecomp_nll)) {
    std::cerr << std::setprecision(17)
              << "FAIL: agecomp_nll got " << likelihood.agecomp_nll
              << " expected " << expected_agecomp_nll
              << " diff " << (likelihood.agecomp_nll - expected_agecomp_nll)
              << "\n";
    return EXIT_FAILURE;
  }

  if (!nearly_equal(likelihood.total_nll, expected_total_nll)) {
    std::cerr << std::setprecision(17)
              << "FAIL: total_nll got " << likelihood.total_nll
              << " expected " << expected_total_nll
              << " diff " << (likelihood.total_nll - expected_total_nll)
              << "\n";
    return EXIT_FAILURE;
  }

  std::cout << "PASSED: Bigeye v2 CAA LikelihoodPackage regression\n";
  return EXIT_SUCCESS;
}
CPP

cat > run_bigeye_v2_14_likelihood_package_caa_check.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples

c++ -std=c++17 -O3 \
  -I. \
  examples/NMFS/pifsc_bigeye_tuna/v2/14_likelihood_package_caa/bigeye_v2_14_likelihood_package_caa_check.cpp \
  -o build/examples/bigeye_v2_14_likelihood_package_caa_check

./build/examples/bigeye_v2_14_likelihood_package_caa_check
SH

chmod +x run_bigeye_v2_14_likelihood_package_caa_check.sh

echo "created CAA LikelihoodPackage"
