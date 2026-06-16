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
