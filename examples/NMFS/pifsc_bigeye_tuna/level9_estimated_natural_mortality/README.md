# Level 9: Estimated Natural Mortality

Level 9 keeps the best-supported Level 6 structure and adds natural mortality
as an estimated fixed effect.

This is a diagnostic experiment. Estimating M is often weakly identified, but
letting it move can reveal whether fixed survival assumptions are forcing
recruitment deviations to compensate for age-structure mismatch.

## Compared with Level 6

- adds `log_m`
- keeps fleet-specific q
- keeps longline logistic selectivity
- keeps purse-seine fixed age-based selectivity
- keeps annual recruitment deviations
- keeps recruitment and safe wiggle diagnostics
