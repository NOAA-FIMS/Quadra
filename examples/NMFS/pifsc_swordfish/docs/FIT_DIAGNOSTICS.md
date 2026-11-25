# Swordfish fit diagnostics and management benchmarks

These checks use the joint penalized likelihood with 84 independent fitting
coordinates. Recruitment deviations are fitted jointly; this is the SS3
comparison target, with fixed biology, steepness, sigmaR, and bias adjustment.
The calculations do not integrate recruitment deviations with Laplace.

## Reproduce

After `make fit-swordfish`, run:

```
make uncertainty
make robustness
make test-benchmarks
```

A reviewable snapshot is saved in [reference/fit_diagnostics](../reference/fit_diagnostics).
Outputs go in `build/`; the canonical fitted parameters are read without being
overwritten. Uncertainty requires Eigen, already provided by Quadra. Python
comparison scripts use only the standard library.

## Curvature and uncertainty

`hessian_coarse.csv` and `hessian_fine.csv` contain central differences of exact
AD gradients at logit-coordinate steps 0.0001 and 0.00005. Rows identify the
perturbed coordinate; columns identify the gradient component. The symmetric
fine-step matrix is inverted only after checking positive definiteness,
relative asymmetry, and agreement between steps. No ridge or eigenvalue
truncation is applied.

The minimum eigenvalue is approximately 0.00194346, the maximum 40642.6,
and the condition number 20.9 million. Relative matrix asymmetry is 1.61e-7,
relative difference between step sizes 2.34e-6, and the largest relative
parameter SE change is 9.49e-6.

`parameter_uncertainty.csv` reports all 84 independent parameter estimates,
delta-method native-scale standard errors, and 95% coordinate-Wald intervals
transformed through each parameter's bounded logit. These intervals are
asymmetric on the native scale. The dependent final main recruitment deviation,
profiled Q, and fixed parameters are not extra independent coordinates.
`coordinate_covariance.csv` preserves the full covariance for derived quantities.

`derived_uncertainty.csv` uses numerical derivatives of derived quantities and
the full covariance. Its 95% intervals are symmetric delta-method intervals.
These are local likelihood approximations conditional on this model setup;
they do not include structural or scenario uncertainty. All 28 estimated
control-parameter SEs and 11 selected derived-quantity SEs are compared with
the supplied report by `compare_uncertainty.py` (1% relative tolerance). All 39
comparisons pass; the maximum relative difference is 9.30e-5 (0.0093%). Halving
the derived-quantity differentiation step changes its SEs by at most 6.84e-6
relative. For example, terminal SSB is 35,778.13 t with SE 6,008.87 t; MSY is
14,924.49 t with SE 648.44 t.

## Perturbed starts

`robustness.csv` records three reproducible fits (seeds 1701–1703), with independent
normal perturbations of every original logit coordinate at SDs 0.02, 0.05,
and 0.10 respectively. These are local robustness checks around the fitted
solution, not a global-optimum proof or starts from unfitted control values.

Optimization uses a constant Cholesky factor of the fitted coordinate covariance
to rescale coordinates. This is a linear change of variables with no change to
the likelihood, bounds, or MLE target; no density Jacobian is added. Each run
checks the objective and AD gradient chain rule in the new coordinates.
Convergence requires both the optimizer's rescaled tolerance and an independently
evaluated original-coordinate gradient norm below 0.0001. The CSV records the
original-coordinate gradient norm and distance from the canonical fit. All three
starts converged: their initial objectives were 300.327, 345.868, and 586.043;
their final objectives agree with the baseline within 1.97e-10. Original-coordinate
gradient norms were 1.60e-6, 5.95e-7, and 1.68e-6.

## Equilibrium reference points

`management_reference_points.csv` is calculated independently of SS3 target
values. It holds the fitted 2021 fleet/season fishing pattern constant and scales
it uniformly. Biology and selectivity use the configured terminal benchmark
year. Beverton–Holt equilibrium recruitment is solved analytically for each
fishing scale, SPR and spawning-biomass targets use bracketed bisection, and MSY
uses a bracketed root of a fourth-order numerical yield derivative. The derivative
root avoids rounding noise near a flat yield maximum when differentiating MSY.

Spawning timing and annual F follow the explicit SS3 v3.30.20 conventions in
`MODEL_AND_COMPARISON.md`. Numbers are thousands of fish, biomass and dead yield
are tonnes, and annual F uses SS3 units 4 pooled across ages 1–10.

| Quantity | Quadra |
|---|---:|
| Unfished SSB | 99,807.76 t |
| SSB at 20% biomass target | 19,961.55 t |
| F at 20% biomass target | 0.158320 |
| Yield at 20% biomass target | 14,814.71 t |
| SSB at SPR 20% | 17,680.23 t |
| F at SPR 20% | 0.170494 |
| Yield at SPR 20% | 14,908.70 t |
| SSB at MSY | 16,387.88 t |
| FMSY | 0.178078 |
| MSY | 14,924.49 t |

All ten comparisons pass against the rounded SS3 report with a relative
threshold of 1e-5. The MSY SSB difference includes the precision of SS3's
benchmark search. Unit checks also exercise additional target values and verify
that yield is lower on both sides of the MSY solution.

Dynamic unfished biomass now matches SS3; deterministic scenarios and their
comparison with published projections are documented in [PROJECTIONS.md](PROJECTIONS.md).
Published bootstrap replay, profile likelihoods, and broader model sensitivities
remain follow-up work. The supplied
forecast configuration requests no projections. The striped marlin model still
requires its own inputs and reference model.
