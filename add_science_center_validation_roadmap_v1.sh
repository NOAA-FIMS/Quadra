#!/usr/bin/env bash
set -euo pipefail

echo "== Add science-center validation roadmap and SEFSC red-snapper scaffold =="

mkdir -p docs/validation
mkdir -p examples/NMFS/sefsc_red_snapper/{data,quadra,tmb,outputs,validation}

cat > docs/validation/science-center-example-roadmap.md <<'MD'
# Science Center Example Validation Roadmap

This document tracks a proposed validation suite with one representative assessment-style example from each NOAA Fisheries Science Center.

The goal is to build examples that are:
- public-data-safe or synthetic,
- reproducible,
- paired with TMB reference implementations where practical,
- documented with expected outputs,
- capable of reporting uncertainty, derived quantities, and projections.

## Proposed example set

| Science Center | Example | Status | Main validation target |
|---|---|---:|---|
| PIFSC | Opakapaka projection example | In progress | Projection validation and Level-1 uncertainty reporting |
| SEFSC | Red-snapper-style age-structured model | Scaffolded | Age structure, selectivity, recruitment deviations, projections |
| NEFSC | Groundfish/index-heavy assessment | Planned | Multiple indices, survey likelihoods, retrospective-style diagnostics |
| NWFSC | West Coast age-structured model | Planned | Age composition, selectivity, biological reference points |
| AFSC | Pollock/sablefish-style model | Planned | Recruitment deviations, state-space/random-effect scalability |
| SWFSC | CPS/tuna-style model | Planned | Time-varying dynamics, index scaling, projection scenarios |

## Shared validation requirements

Each example should eventually include:

1. Quadra implementation
2. TMB comparison implementation
3. synthetic or public-data-safe input data
4. reproducible runner
5. fit diagnostics
6. standard errors and confidence intervals
7. random-effect conditional uncertainty
8. derived quantity uncertainty
9. projection envelopes
10. comparison summary against TMB

## Recommended directory layout

```text
examples/<example_name>/
  README.md
  data/
  quadra/
  tmb/
  outputs/
  validation/
```

## Development order

1. Finish Opakapaka Level-1 uncertainty reporting.
2. Scaffold SEFSC red-snapper-style age-structured model.
3. Add minimal Quadra implementation.
4. Add TMB reference implementation.
5. Add validation summary and uncertainty outputs.
6. Repeat for the remaining science centers.
MD

cat > examples/NMFS/sefsc_red_snapper/README.md <<'MD'
# SEFSC Red-Snapper-Style Assessment Example

This directory is a placeholder for a synthetic, public-data-safe red-snapper-style assessment example.

The goal is not to reproduce an official assessment. The goal is to provide a representative SEFSC-style validation case for Quadra with age structure, selectivity, recruitment deviations, uncertainty reporting, and projections.

## Planned model features

- age-structured population dynamics
- catch likelihood
- survey/index likelihood
- age-composition likelihood
- recruitment deviations as random effects
- age-based selectivity
- derived quantities:
  - biomass
  - spawning biomass proxy
  - depletion
  - fishing mortality proxy
  - MSY-like reference metrics
- projection scenarios
- uncertainty outputs:
  - inverse Hessian / covariance
  - standard errors
  - confidence intervals
  - random-effect conditional uncertainty
  - derived quantity uncertainty
  - projection envelopes

## Directory layout

```text
data/        synthetic or public-data-safe inputs
quadra/      Quadra implementation
tmb/         TMB reference implementation
outputs/     generated outputs, ignored by git
validation/  comparison summaries and validation notes
```

## Initial validation target

The first milestone is a minimal working model with:

1. deterministic age-structured dynamics,
2. one abundance index,
3. synthetic catch observations,
4. recruitment deviations,
5. TMB side-by-side comparison,
6. Level-1 uncertainty outputs.
MD

cat > examples/NMFS/sefsc_red_snapper/validation/validation_plan.md <<'MD'
# SEFSC Red-Snapper-Style Validation Plan

## Level 0: deterministic fit

- Build a minimal deterministic age-structured model.
- Fit fixed effects only.
- Confirm objective value and parameter estimates are stable.

## Level 1: random effects and uncertainty

- Add recruitment deviations as random effects.
- Extract conditional random-effect uncertainty.
- Add fixed-effect covariance and confidence intervals.
- Add derived quantity uncertainty.

## Level 2: TMB comparison

- Implement matching TMB reference model.
- Compare:
  - objective value
  - fixed-effect estimates
  - random-effect modes
  - standard errors
  - derived quantities
  - projection summaries

## Level 3: projections

- Add projection scenarios.
- Report projection envelopes.
- Compare Quadra and TMB projection outputs where feasible.

## Notes

This example should remain synthetic or public-data-safe. It should not be presented as an official red snapper assessment.
MD

cat > examples/NMFS/sefsc_red_snapper/data/README.md <<'MD'
# Data

Synthetic or public-data-safe input files will live here.

Do not commit generated outputs or confidential assessment data.
MD

cat > examples/NMFS/sefsc_red_snapper/quadra/README.md <<'MD'
# Quadra Implementation

Quadra model source files for the SEFSC red-snapper-style example will live here.
MD

cat > examples/NMFS/sefsc_red_snapper/tmb/README.md <<'MD'
# TMB Reference Implementation

TMB comparison files for the SEFSC red-snapper-style example will live here.
MD

cat > examples/NMFS/sefsc_red_snapper/outputs/.gitignore <<'EOF'
*
!.gitignore
EOF

echo
echo "Created:"
echo "  docs/validation/science-center-example-roadmap.md"
echo "  examples/NMFS/sefsc_red_snapper/"
echo
echo "Next:"
echo "  git add docs/validation/science-center-example-roadmap.md examples/NMFS/sefsc_red_snapper"
echo "  git commit -m \"Add science center validation roadmap and SEFSC scaffold\""
