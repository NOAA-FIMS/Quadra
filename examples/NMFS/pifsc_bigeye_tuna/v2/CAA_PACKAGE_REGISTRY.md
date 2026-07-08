# CAA Package Registry

Generated from `architecture/packages/*/package.meta`.

| Package | Consumes | Produces | Steps |
|---|---|---|---|
| LifeHistoryPackage | BigeyeModelData, LifeHistoryParameters | LifeHistoryState | BigeyeLifeHistory |
| PopulationPackage | BigeyeModelData, PopulationParameters, LifeHistoryState, FleetState | PopulationState | FixedRecruitment<br>Survival<br>Aging<br>PlusGroup<br>SpawningBiomass |
| MovementPackage | BigeyeModelData, MovementParameters, PopulationState | PopulationState | IdentityMovement |
| FleetPackage | BigeyeModelData, FleetParameters, LifeHistoryState, PopulationState | FleetState | LogisticSelectivity<br>FishingMortality<br>BaranovCatch |
| ObservationPackage | BigeyeModelData, FleetParameters, PopulationState, FleetState | FleetState | BiomassIndexPrediction<br>CatchAgeCompositionPrediction |
| LikelihoodPackage | BigeyeModelData, FleetParameters, FleetState | LikelihoodState | LognormalCatchLikelihood<br>LognormalIndexLikelihood<br>MultinomialAgeCompLikelihood |
