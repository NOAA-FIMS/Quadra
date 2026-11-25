# Swordfish projection parameter uncertainty

All five 2022–2031 scenarios now have local parameter uncertainty for spawning
biomass, recruitment, January biomass, dead catch and annual F. The 250 point
estimates reproduce the deterministic forecasts. Numerical validation checks
propagation through the complete scenario calculation.

[Five-page chart](../reference/projection_uncertainty/projection_uncertainty.pdf)
(one page per quantity), [all intervals](../reference/projection_uncertainty/projection_uncertainty.csv),
and [2031 results](../reference/projection_uncertainty/projection_uncertainty_2031.csv).

## Reproduce

From this example directory, after fitting and preparing parameter uncertainty:

```
make projection-uncertainty
make plot-projection-uncertainty
```

If the parameter covariance does not exist or predates the fitted parameter CSV,
the first target generates it using the existing `uncertainty` workflow. Covariance provenance is now mandatory; artifacts from earlier workflow versions
without a stamp must be regenerated. To
force a new curvature calculation after changing the model code, run
`make uncertainty` first. Plotting additionally requires base R.

The numerical unit check is available separately:

```
make build/test_delta_method
./build/test_delta_method
```

`make test-projection-uncertainty` runs that check together with the full
projection calculation and output validation. Outputs are in `build/projections/`;
a validated snapshot is in `reference/projection_uncertainty/`.

## What varies

The input is the full covariance of 84 independent bounded-logit coordinates
from the joint penalized likelihood. It includes correlations among estimated
control parameters and historical recruitment deviations. The final constrained
main recruitment deviation is derived by the same zero-sum rule as in fitting;
it is not given an independent variance.

Each parameter perturbation recomputes the historical population, January 2022
starting state, dynamic unfished baseline, 2019–2021 fleet/season fishing mixture,
and all five scenario targets. Historical SPR targets retain the historical
selectivity blocks. The complete forecast trajectory is then recalculated.
Deterministic and uncertainty calculations share the scenario parser and resolver
in `include/swordfish/scenarios.hpp`.

This treats each scenario as a **rule recalibrated for a possible parameter
value**. It does not hold numerical fleet fishing rates fixed at their fitted
values. Fixed-rate uncertainty would answer a different policy question. Target
recalibration can offset some parameter effects and narrow intervals, particularly
for catches constrained by historical fishing/catch relationships.

Biology and hyperparameters fixed in the fitted model remain fixed. Recruitment
after 2021 follows deterministic Beverton–Holt expectation with no future
deviations or log-bias adjustment. No future recruitment process variation,
structural model uncertainty, or alternative scenario descriptions are added.
The separate interpretation sensitivities remain in `PROJECTIONS.md`.

## Calculation and interval interpretation

For the vector of projected outputs `g(x)`, central differences at logit steps
0.0001 and 0.00005 produce two Jacobians. The fine-step Jacobian propagates the
full covariance as `J V J'`. Output covariance retains correlations across
quantities, years and scenarios; parameter correlations are not discarded.

Native-scale standard errors are the square roots of its diagonal. To keep
intervals positive, the pointwise 95% interval for estimate `g` is
`g * exp(+-1.959964 * SE/g)`: a first-order log-scale Wald approximation. These
are neither simultaneous trajectory bands nor intervals for realized future
outcomes. They are not the assessment's bootstrap medians or bootstrap intervals,
and do not establish probabilities of crossing management thresholds.

The local covariance approximation has not been calibrated for frequentist
coverage by repeated simulations. The saved checks establish numerical stability
and consistency, not exact coverage. Existing differences from published Table 12
remain unresolved; adding parameter intervals does not reproduce its bootstrap
analysis.

## Selected 2031 results

All figures below are tonnes; intervals are pointwise local 95% intervals under
the recalibrated target rules described above.

| Scenario | SSB | SSB interval | Dead catch | Catch interval |
|---|---:|---:|---:|---:|
| Dynamic20 | 19,446 | 17,130–22,075 | 15,083 | 13,674–16,637 |
| Reference2008_2010 | 31,470 | 22,566–43,889 | 13,362 | 12,706–14,051 |
| SPR30 | 27,674 | 24,543–31,203 | 14,031 | 12,713–15,485 |
| MSY | 16,780 | 14,886–18,916 | 15,247 | 13,817–16,825 |
| StatusQuo2019_2021 | 38,063 | 27,961–51,815 | 12,018 | 11,465–12,598 |

These differences in interval width are conditional on each target rule. They
should not be interpreted as different levels of total ecological uncertainty.
Recruitment, biomass and F intervals are included in the full CSV and chart.

## Validation and saved artifacts

- All 250 point estimates agree with deterministic projections within 7.11e-15
  relative. Refactoring the common scenario code left the saved deterministic
  forecast outputs unchanged.
- Halving the differentiation step changes standard errors by at most 4.56e-6
  relative (0.000456%). The largest relative Jacobian row change is 5.55e-7.
- Four independent directional differences agree with the composed Jacobian;
  maximum scaled directional error is 1.16e-5. Directional scaling uses the
  larger of the derivative magnitude and 0.001 times the output estimate.
- Scalar recomputations check selected cross-year/cross-quantity covariance
  entries. All output variances equal the reported SEs squared; covariance
  symmetry and all positive interval endpoints are checked.
- A separate analytic test verifies correlated covariance propagation, invariance
  under linear changes of parameter coordinates, zero-variance limits and
  rejection of indefinite parameter covariance.
- All 1,551 dynamic unfished and 47 historical SPR comparisons still pass. The
  projection equilibrium/stationarity checks and five Python input tests pass.
  The eight supplied source file hashes remain unchanged.

The reader checks covariance dimensions, coordinate order, positive definiteness
and symmetry, and agreement of the fitted coordinate vector and bounds with the
curvature inputs. Validation also verifies the original source hashes. Input
SHA-256 hashes are recorded in `projection_uncertainty_validation.json`.

The snapshot contains all intervals, numerical diagnostics, the 250-by-84
Jacobian, the 250-by-250 native output covariance, and the input coordinate
covariance and fitted values. The output covariance can be singular because 250
derived quantities depend on only 84 independent coordinates; it is not inverted.
