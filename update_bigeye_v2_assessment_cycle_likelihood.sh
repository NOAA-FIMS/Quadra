#!/usr/bin/env bash
set -euo pipefail

p="examples/NMFS/pifsc_bigeye_tuna/v2/architecture/assessment/assessment_cycle.hpp"

python3 - <<'PY'
from pathlib import Path
p = Path("examples/NMFS/pifsc_bigeye_tuna/v2/architecture/assessment/assessment_cycle.hpp")
s = p.read_text()

s = s.replace(
'#include "../packages/life_history/life_history_package.hpp"\n',
'#include "../packages/life_history/life_history_package.hpp"\n'
'#include "../packages/likelihood/likelihood_context.hpp"\n'
'#include "../packages/likelihood/likelihood_package.hpp"\n'
)

s = s.replace(
'''    ObservationContext<T> observation_context{
        &parameters.fleets[0],
        &state.populations[0],
        &state.fleets[0]};''',
'''    ObservationContext<T> observation_context{
        &parameters.fleets[0],
        &state.populations[0],
        &state.fleets[0]};

    LikelihoodContext<T> likelihood_context{
        &parameters.fleets[0],
        &state.fleets[0],
        &state.likelihood};'''
)

s = s.replace(
'''    // Prediction-only observation layer.
    ObservationPackage{}(data, observation_context);''',
'''    // Prediction layer.
    ObservationPackage{}(data, observation_context);

    // Objective layer.
    LikelihoodPackage{}(data, likelihood_context);'''
)

p.write_text(s)
PY

echo "updated AssessmentCycle to include LikelihoodPackage"
