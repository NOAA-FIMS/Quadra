
# Quadra Testing Plan

## Unit-style model tests

### CurvatureDependsOnTheta

Purpose:
- validates Laplace log-det gradient contribution,
- validates exact `Hdot`,
- has an analytic optimum.

Expected optimum:

```text
theta* = -0.5 * n_random
```

### PoissonRandomIntercept

Purpose:
- tests nonlinear likelihood,
- tests iid random effects,
- exercises non-Gaussian observation model.

### GaussianRandomWalk

Purpose:
- tests banded sparse Hessian structure,
- tests temporal random effects,
- resembles state-space/random-walk fisheries models.

## Hdot validation

Compile with:

```bash
-DQUADRA_VALIDATE_HDOT
```

Expected for curvature-dependent models:

```text
exact_norm ~= fd_norm
rel_err small
```

Do not run validation on large random-effect systems.

## Fisheries examples

### fisheries_random_year_effects

CPUE/index-style observation with latent year effects.

### fisheries_age_selectivity_random_walk

Logistic selectivity curve with random-walk deviations over age.

### fisheries_index_cpue_laplace

Latent abundance process model with CPUE observations.

## Next tests to add

- Dirichlet-multinomial likelihood
- multinomial age composition
- Beverton-Holt recruitment random effects
- AR1 process prior
- separable fleet/year likelihood
- multi-fleet index model
- surplus production model with process error
