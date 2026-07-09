# CAA Operation Catalog

Generated from `architecture/packages/*/package.meta`.

## LifeHistoryPackage

**Purpose:** Computes life-history quantities

**Operations:**

- ComputeNaturalMortality
- ComputeWeightAtAge
- ComputeMaturityAtAge

---

## PopulationPackage

**Purpose:** Composes population dynamics for one population

**Operations:**

- Recruitment
- Survival
- Aging
- PlusGroup
- SpawningBiomass

---

## MovementPackage

**Purpose:** Composes movement across populations

**Operations:**

- IdentityMovement

---

## FleetPackage

**Purpose:** Composes fleet dynamics for one fleet

**Operations:**

- LogisticSelectivity
- FishingMortality
- BaranovCatch

---

## ObservationPackage

**Purpose:** Composes predicted observations for one fleet

**Operations:**

- BiomassIndexPrediction
- CatchAgeCompositionPrediction

---

## LikelihoodPackage

**Purpose:** Composes likelihood components for one fleet

**Operations:**

- LognormalCatchLikelihood
- LognormalIndexLikelihood
- MultinomialAgeCompLikelihood

---
