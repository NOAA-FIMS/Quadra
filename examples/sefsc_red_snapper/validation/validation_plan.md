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
