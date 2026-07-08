#!/usr/bin/env bash
set -euo pipefail

BASE="examples/NMFS/pifsc_bigeye_tuna/v2/architecture/packages"

cat > "$BASE/life_history/README.md" <<'MD'
# Life History Package

Computes life-history quantities.

## State

- `LifeHistoryState`

## Parameters

- `LifeHistoryParameters`

## Steps

- `BigeyeLifeHistory`

## Produces

- natural mortality-at-age
- weight-at-age
- maturity-at-age
MD

cat > "$BASE/population/README.md" <<'MD'
# Population Package

Composes population dynamics for one population.

## State

- `PopulationState`

## Parameters

- `PopulationParameters`

## Steps

- `FixedRecruitment`
- `Survival`
- `Aging`
- `PlusGroup`
- `SpawningBiomass`

## Produces

- recruits by year
- survivors-at-age
- numbers-at-age
- spawning biomass by year
MD

cat > "$BASE/fleet/README.md" <<'MD'
# Fleet Package

Composes fleet dynamics for one fleet.

## State

- `FleetState`

## Parameters

- `FleetParameters`

## Steps

- `LogisticSelectivity`
- `FishingMortality`
- `BaranovCatch`

## Produces

- selectivity-at-age
- fishing mortality-at-age
- total mortality-at-age
- catch numbers-at-age
- catch biomass
MD

cat > "$BASE/movement/README.md" <<'MD'
# Movement Package

Composes movement across populations.

## State

- `PopulationState`

## Parameters

- `MovementParameters`

## Steps

- `IdentityMovement`

## Produces

- updated population numbers-at-age

Current implementation is identity movement.
MD

cat > "$BASE/observation/README.md" <<'MD'
# Observation Package

Composes predicted observations for one fleet.

## State

- `FleetState`
- `PopulationState`

## Parameters

- `FleetParameters`

## Steps

- `BiomassIndexPrediction`
- `CatchAgeCompositionPrediction`

## Produces

- predicted biomass index
- predicted catch age composition
MD

cat > "$BASE/likelihood/README.md" <<'MD'
# Likelihood Package

Composes likelihood components for one fleet.

## State

- `FleetState`
- `LikelihoodState`

## Parameters

- `FleetParameters`

## Steps

- `LognormalCatchLikelihood`
- `LognormalIndexLikelihood`
- `MultinomialAgeCompLikelihood`

## Produces

- catch negative log-likelihood
- index negative log-likelihood
- age composition negative log-likelihood
- total negative log-likelihood
MD

echo "added CAA package READMEs"
