# TMB Random-Intercept Comparison

This benchmark is intended to compare Quadra and TMB on the same random-intercept likelihood.

Initial scaffold:

- `quadra_random_intercept_compare.cpp` — Quadra implementation.
- `random_intercept_tmb.cpp` — TMB objective template.
- `run_tmb_random_intercept.R` — TMB runner.
- `comparison_outputs/` — CSV outputs.

The TMB benchmark is optional and requires an R environment with TMB installed.

