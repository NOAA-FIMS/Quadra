#!/usr/bin/env bash
set -euo pipefail

BASE="examples/NMFS/pifsc_bigeye_tuna/v2"

cat > "$BASE/architecture/packages/fleet/fleet_context.hpp" <<'CPP'
#pragma once

#include "../../parameters/fleet_parameters.hpp"
#include "../../state/fleet_state.hpp"
#include "../../state/life_history_state.hpp"
#include "../../state/population_state.hpp"

namespace bigeye_v2 {

template <typename T>
struct FleetContext {
  const FleetParameters<T> *parameters = nullptr;
  const LifeHistoryState<T> *life = nullptr;
  const PopulationState<T> *population = nullptr;
  FleetState<T> *fleet = nullptr;
};

} // namespace bigeye_v2
CPP

python3 - <<'PY'
from pathlib import Path

p = Path("examples/NMFS/pifsc_bigeye_tuna/v2/architecture/packages/fleet/fleet_package.hpp")
s = p.read_text()

if '#include "fleet_context.hpp"' not in s:
    s = s.replace('#pragma once\n\n', '#pragma once\n\n#include "fleet_context.hpp"\n')

old = '''  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const FleetParameters<T> &parameters,
                  const LifeHistoryState<T> &life,
                  const PopulationState<T> &population,
                  FleetState<T> &fleet) const {
    LogisticSelectivity{}(data, parameters, fleet);
    FishingMortality{}(parameters, life, fleet);
    BaranovCatch{}(data, life, population, fleet);
  }'''

new = old + '''

  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const FleetContext<T> &context) const {
    LogisticSelectivity{}(data, *context.parameters, *context.fleet);
    FishingMortality{}(*context.parameters, *context.life, *context.fleet);
    BaranovCatch{}(data, *context.life, *context.population, *context.fleet);
  }'''

if 'FleetContext<T> &context' not in s:
    s = s.replace(old, new)

p.write_text(s)
PY

cat > "$BASE/architecture/packages/population/population_context.hpp" <<'CPP'
#pragma once

#include "../../parameters/population_parameters.hpp"
#include "../../state/fleet_state.hpp"
#include "../../state/life_history_state.hpp"
#include "../../state/population_state.hpp"

namespace bigeye_v2 {

template <typename T>
struct PopulationContext {
  const PopulationParameters<T> *parameters = nullptr;
  const LifeHistoryState<T> *life = nullptr;
  const FleetState<T> *fleet = nullptr;
  PopulationState<T> *population = nullptr;
};

} // namespace bigeye_v2
CPP

python3 - <<'PY'
from pathlib import Path

p = Path("examples/NMFS/pifsc_bigeye_tuna/v2/architecture/packages/population/population_package.hpp")
s = p.read_text()

if '#include "population_context.hpp"' not in s:
    s = s.replace('#pragma once\n\n', '#pragma once\n\n#include "population_context.hpp"\n')

old = '''  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const PopulationParameters<T> &parameters,
                  const LifeHistoryState<T> &life,
                  const FleetState<T> &fleet,
                  PopulationState<T> &population) const {'''

if 'PopulationContext<T> &context' not in s:
    insert = '''
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const PopulationContext<T> &context) const {
    (*this)(data,
            *context.parameters,
            *context.life,
            *context.fleet,
            *context.population);
  }

'''
    s = s.replace('struct PopulationPackage {\n', 'struct PopulationPackage {\n' + insert)

p.write_text(s)
PY

cat > "$BASE/architecture/packages/observation/observation_context.hpp" <<'CPP'
#pragma once

#include "../../parameters/fleet_parameters.hpp"
#include "../../state/fleet_state.hpp"
#include "../../state/population_state.hpp"

namespace bigeye_v2 {

template <typename T>
struct ObservationContext {
  const FleetParameters<T> *parameters = nullptr;
  const PopulationState<T> *population = nullptr;
  FleetState<T> *fleet = nullptr;
};

} // namespace bigeye_v2
CPP

python3 - <<'PY'
from pathlib import Path

p = Path("examples/NMFS/pifsc_bigeye_tuna/v2/architecture/packages/observation/observation_package.hpp")
s = p.read_text()

if '#include "observation_context.hpp"' not in s:
    s = s.replace('#pragma once\n\n', '#pragma once\n\n#include "observation_context.hpp"\n')

if 'ObservationContext<T> &context' not in s:
    insert = '''
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const ObservationContext<T> &context) const {
    BiomassIndexPrediction{}(data,
                             *context.parameters,
                             *context.population,
                             *context.fleet);
    CatchAgeCompositionPrediction{}(data, *context.fleet);
  }

'''
    s = s.replace('struct ObservationPackage {\n', 'struct ObservationPackage {\n' + insert)

p.write_text(s)
PY

cat > "$BASE/architecture/packages/likelihood/likelihood_context.hpp" <<'CPP'
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
CPP

python3 - <<'PY'
from pathlib import Path

p = Path("examples/NMFS/pifsc_bigeye_tuna/v2/architecture/packages/likelihood/likelihood_package.hpp")
s = p.read_text()

if '#include "likelihood_context.hpp"' not in s:
    s = s.replace('#pragma once\n\n', '#pragma once\n\n#include "likelihood_context.hpp"\n')

if 'LikelihoodContext<T> &context' not in s:
    insert = '''
  template <typename T>
  void operator()(const BigeyeModelData<T> &data,
                  const LikelihoodContext<T> &context) const {
    (*this)(data,
            *context.parameters,
            *context.fleet,
            *context.likelihood);
  }

'''
    s = s.replace('struct LikelihoodPackage {\n', 'struct LikelihoodPackage {\n' + insert)

p.write_text(s)
PY

cat > run_bigeye_v2_package_interface_unification_check.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

./run_bigeye_v2_08_population_package_orchestrator_caa_check.sh
./run_bigeye_v2_09_fleet_package_orchestrator_caa_check.sh
./run_bigeye_v2_13_observation_package_caa_check.sh
./run_bigeye_v2_14_likelihood_package_caa_check.sh
./run_bigeye_v2_10_assessment_cycle_caa_check.sh

echo
echo "PASSED: Bigeye v2 CAA package interface unification check"
SH

chmod +x run_bigeye_v2_package_interface_unification_check.sh

echo "unified package interfaces with context adapters"
