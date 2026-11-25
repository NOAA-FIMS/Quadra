# Reference-point uncertainty in the native model

Run from the swordfish case directory:

```bash
make pifsc_swordfish
./build/pifsc_swordfish --input build/native/swordfish-ss-par-start.json \
  --output build/native/my-reference-uncertainty --uncertainty
```

The model first fits and checks its final gradient, then calculates fresh local
parameter covariance and propagates it to the reference points. The optional
flag keeps ordinary fits fast: the extra work evaluates two finite-difference
steps for all 84 independent coordinates and prints progress for each step.
It requires fit mode and does not read SS3 files, saved covariance or reference
outputs. A failed curvature/stability check fails the requested run.

## Outputs

`management_reference_point_uncertainty.csv` contains one row per benchmark
and quantity: estimate, standard error, CV, lower/upper 95% limits, units and
interval method. It includes all six standard quantities for unfished, Btgt20,
SPR20 and MSY: 24 rows in total.

Positive quantities use pointwise log-scale Wald intervals:
`estimate * exp(+-1.95996398454005 * SE / estimate)`.
Nonconstant SPR proportions use logit-scale intervals confined to [0,1].
Known unfished zeros, unfished SPR=1 and target SPR=.2 have zero SE and
zero-width intervals. CV is blank for zero estimates. Btgt20 raw SPR is evaluated
as `.2/h` in the supported positive Beverton-Holt equilibrium, so a fixed h
does not acquire spurious root-rounding uncertainty.

Additional audit files are:

- `reference_coordinate_covariance.csv`: full independent-coordinate covariance.
- `management_reference_point_covariance.csv`: full derived-output covariance.
- `reference_hessian_coarse.csv` and `reference_hessian_fine.csv`.
- `reference_jacobian_coarse.csv` and `reference_jacobian_fine.csv`.
- `reference_uncertainty_summary.json`: numerical diagnostics and acceptance status.

The native runner also preserves its input JSON, binary hash, parameter estimates,
ordinary reference-point table and run status. Thus covariance is tied to this
run rather than borrowed from an older fit. Derived covariance has deterministic
and dependent rows and can be singular; it is not a matrix to invert.

## Method and validation

The covariance is the inverse Hessian of the fitted joint penalized likelihood
in independent bounded-logit coordinates. Hessians use centered differences of
exact Quadra AD gradients at steps `1e-4` and `5e-5`. Both must be finite,
symmetric/stable to `1e-5` relative norm and positive definite. Parameter SEs
must agree within 0.1%. No ridge, clipped eigenvalues or diagonal-only covariance
is substituted for an invalid fit.

Reference-point Jacobians are centered differences through the actual equilibrium
root/MSY calculations at both steps. Propagation uses `J V J-transpose`, retaining
all estimated-parameter correlations and the recruitment sum constraint. Reference
SEs must change by less than 0.1% when varying either the Hessian or Jacobian
step. These are numerical derivatives through reference-point solves, not exact
AD derivatives of the roots.

These are conditional local parameter-confidence intervals. Fixed biology and
other fixed settings remain fixed. Historical recruitment deviations are jointly
estimated, not Laplace-integrated. No future process variability, model-selection
uncertainty, bootstrap or simultaneous coverage is added. The implementation uses
the existing correlated delta-method utilities and checks against the earlier
canonical-fit SEs; it does not claim equivalence to a different FIMS likelihood.

The [saved native run](../reference/reference_point_uncertainty/reference_uncertainty_summary.json)
passes both step checks. Its ten corresponding SEs also match the saved printed
SS3 SEs, with [maximum relative difference below 0.0006%](../reference/reference_point_uncertainty/SS3_reference_SE_comparison.json).
For example:

| MSY quantity | Estimate | SE | CV | 95% interval |
| --- | ---: | ---: | ---: | --- |
| Dead catch, tonnes/year | 14,924.50 | 648.44 | 4.34% | 13,706.18-16,251.10 |
| Female SSB, tonnes | 16,387.89 | 813.98 | 4.97% | 14,867.71-18,063.51 |
| Annual F | 0.178078 | 0.003997 | 2.24% | 0.170414-0.186086 |

These bounds use the log-scale method above; earlier canonical diagnostics used
symmetric native-scale Wald intervals, so SEs agree while interval endpoints
can differ. The [PDF figure](../output/pdf/swordfish_reference_point_uncertainty.pdf)
shows the positive pointwise intervals.

For an optional scientific PDF (requires base R):

```bash
Rscript scripts/plot_reference_uncertainty.R \
  build/native/my-reference-uncertainty/management_reference_point_uncertainty.csv \
  output/pdf/swordfish_reference_point_uncertainty.pdf
```
