#!/usr/bin/env bash
set -euo pipefail

BASE="examples/NMFS/pifsc_bigeye_tuna/v2"

cat > "$BASE/caa_manifest.yml" <<'YAML'
assessment:
  name: Bigeye Tuna Reference Assessment
  architecture: Composable Assessment Architecture
  version: 0.2

cycle:
  - LifeHistoryPackage
  - PopulationPackage
  - MovementPackage
  - FleetPackage
  - ObservationPackage
  - LikelihoodPackage

packages:
  life_history:
    package: LifeHistoryPackage
    step: BigeyeLifeHistory

  population:
    package: PopulationPackage
    steps:
      - FixedRecruitment
      - Survival
      - Aging
      - PlusGroup
      - SpawningBiomass

  movement:
    package: MovementPackage
    step: IdentityMovement

  fleet:
    package: FleetPackage
    steps:
      - LogisticSelectivity
      - FishingMortality
      - BaranovCatch

  observation:
    package: ObservationPackage
    steps:
      - BiomassIndexPrediction
      - CatchAgeCompositionPrediction

  likelihood:
    package: LikelihoodPackage
    steps:
      - LognormalCatchLikelihood
      - LognormalIndexLikelihood
      - MultinomialAgeCompLikelihood
YAML

cat > inspect_bigeye_v2_caa_manifest.sh <<'SH'
#!/usr/bin/env bash
set -euo pipefail

MANIFEST="examples/NMFS/pifsc_bigeye_tuna/v2/caa_manifest.yml"

echo "CAA Manifest"
echo
grep -E "name:|architecture:|version:" "$MANIFEST"
echo
echo "Assessment cycle:"
awk '
  /^cycle:/ {in_cycle=1; next}
  /^packages:/ {in_cycle=0}
  in_cycle && /^  - / {print "  " $0}
' "$MANIFEST"
SH

chmod +x inspect_bigeye_v2_caa_manifest.sh

echo "created CAA v0.2 manifest"
