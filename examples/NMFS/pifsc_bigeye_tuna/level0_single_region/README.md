# Level 0: Single-Region Bigeye Tuna Prototype

This is the baseline model in the diagnostic-guided modeling ladder.

## Purpose

Establish a minimal age-structured tuna-like model before adding fleets,
regions, movement, or tagging.

## Intended structure

- one region
- one aggregate catch series
- one abundance index
- age-structured population dynamics
- recruitment deviations as random effects
- fixed effects for recruitment scale, fishing mortality, and catchability

## Diagnostic questions

- Does the model converge cleanly?
- Is the random-effect Hessian positive definite?
- Are recruitment deviations weakly or strongly correlated?
- Is curvature concentrated in a few years?
- Are catchability, recruitment scale, and fishing mortality separable?
- What diagnostics worsen when the next level adds fleets?

## Status

Scaffold only. The first implementation should borrow the Red Snapper layout,
then rename/reparameterize toward a tuna-like example.
