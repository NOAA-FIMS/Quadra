# Validated local projection uncertainty snapshot

Methods, interpretation, and reproduction commands:
[PROJECTION_UNCERTAINTY.md](../../PROJECTION_UNCERTAINTY.md).

Contains 250 pointwise intervals across five scenarios, ten years and five
quantities, with the complete Jacobian and native output covariance. Input
coordinate covariance, coordinate mapping and fitted values are retained.
`projection_uncertainty_validation.json` records input hashes and numerical
checks. All five PDF chart pages were rendered and visually reviewed.

Scenario targets and fleet mixtures are recalibrated under each parameter
perturbation. Intervals reflect local parameter uncertainty only. They exclude
future recruitment process variation and do not replicate the assessment's
bootstrap forecast results.
