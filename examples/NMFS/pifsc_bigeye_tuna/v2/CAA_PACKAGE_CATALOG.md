# CAA Package Catalog

Generated from `architecture/packages/*/package.meta`.

## LifeHistoryPackage

**Purpose:** Computes life-history quantities

**Consumes:** BigeyeModelData, LifeHistoryParameters

**Produces:** LifeHistoryState

**Steps:**

- BigeyeLifeHistory

---

## PopulationPackage

**Purpose:** Composes population dynamics for one population

**Consumes:** BigeyeModelData, PopulationParameters, LifeHistoryState, FleetState

**Produces:** PopulationState

**Steps:**

- FixedRecruitment
- Survival
- Aging
- PlusGroup
- SpawningBiomass

---

## MovementPackage

**Purpose:** Composes movement across populations

**Consumes:** BigeyeModelData, MovementParameters, PopulationState

**Produces:** PopulationState

**Steps:**

- IdentityMovement

---

## FleetPackage

**Purpose:** Composes fleet dynamics for one fleet

**Consumes:** BigeyeModelData, FleetParameters, LifeHistoryState, PopulationState

**Produces:** FleetState

**Steps:**

- LogisticSelectivity
- FishingMortality
- BaranovCatch

---

## ObservationPackage

**Purpose:** Composes predicted observations for one fleet

**Consumes:** BigeyeModelData, FleetParameters, PopulationState, FleetState

**Produces:** FleetState

**Steps:**

- BiomassIndexPrediction
- CatchAgeCompositionPrediction

---

## LikelihoodPackage

**Purpose:** Composes likelihood components for one fleet

**Consumes:** BigeyeModelData, FleetParameters, FleetState

**Produces:** LikelihoodState

**Steps:**

- LognormalCatchLikelihood
- LognormalIndexLikelihood
- MultinomialAgeCompLikelihood

---
