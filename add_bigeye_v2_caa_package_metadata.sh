#!/usr/bin/env bash
set -euo pipefail

BASE="examples/NMFS/pifsc_bigeye_tuna/v2/architecture/packages"

cat > "$BASE/life_history/package.meta" <<'META'
name: LifeHistoryPackage
purpose: Computes life-history quantities
consumes: BigeyeModelData, LifeHistoryParameters
produces: LifeHistoryState
steps:
  - BigeyeLifeHistory
META

cat > "$BASE/population/package.meta" <<'META'
name: PopulationPackage
purpose: Composes population dynamics for one population
consumes: BigeyeModelData, PopulationParameters, LifeHistoryState, FleetState
produces: PopulationState
steps:
  - FixedRecruitment
  - Survival
  - Aging
  - PlusGroup
  - SpawningBiomass
META

cat > "$BASE/movement/package.meta" <<'META'
name: MovementPackage
purpose: Composes movement across populations
consumes: BigeyeModelData, MovementParameters, PopulationState
produces: PopulationState
steps:
  - IdentityMovement
META

cat > "$BASE/fleet/package.meta" <<'META'
name: FleetPackage
purpose: Composes fleet dynamics for one fleet
consumes: BigeyeModelData, FleetParameters, LifeHistoryState, PopulationState
produces: FleetState
steps:
  - LogisticSelectivity
  - FishingMortality
  - BaranovCatch
META

cat > "$BASE/observation/package.meta" <<'META'
name: ObservationPackage
purpose: Composes predicted observations for one fleet
consumes: BigeyeModelData, FleetParameters, PopulationState, FleetState
produces: FleetState
steps:
  - BiomassIndexPrediction
  - CatchAgeCompositionPrediction
META

cat > "$BASE/likelihood/package.meta" <<'META'
name: LikelihoodPackage
purpose: Composes likelihood components for one fleet
consumes: BigeyeModelData, FleetParameters, FleetState
produces: LikelihoodState
steps:
  - LognormalCatchLikelihood
  - LognormalIndexLikelihood
  - MultinomialAgeCompLikelihood
META

echo "added CAA package metadata"
