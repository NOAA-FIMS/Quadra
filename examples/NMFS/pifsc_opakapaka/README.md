# Synthetic opakapaka-style fit and projection example

This example is synthetic and public-data-safe. It is not an official assessment
and does not use confidential or official Pacific Islands data.

The purpose is to demonstrate the intended public Quadra workflow:

```cpp
OpakapakaProjectionModel model(data);

auto params = model.initial_parameters();

auto fit = quadra::optimize_lbfgs(model, params, opts);

auto projection = model.project(fit, projection_options);
```

The example intentionally keeps inference machinery inside Quadra:

- no example-local optimizer
- no finite-difference gradient or Hessian code
- no visible Hessian extraction
- no manual backend selection

`optimize_lbfgs()` returns an `OptResult` containing fixed-effect estimates,
random-effect modes, convergence diagnostics, and structure/backend metadata.

## Run

```bash
./run_opakapaka_projection.sh
```

Outputs are written to:

```text
examples/NMFS/pifsc_opakapaka/outputs/synthetic_fit_summary.csv
examples/NMFS/pifsc_opakapaka/outputs/synthetic_projection_scenarios.csv
```

## Opakapaka Projection Validation

A synthetic opakapaka-style state-space projection model was implemented in both Quadra and TMB using identical data and model structure.

### Numerical Agreement

| Metric | Quadra | TMB |
|----------|----------|----------|
| Objective | 40.774049 | 40.774049 |
| log(q) | -6.505367 | -6.505366 |
| q | 0.001495 | 0.001495 |
| Random effects | 20 | 20 |

The two implementations converge to effectively identical parameter estimates and objective values.

### Quadra Structure Detection

Quadra automatically identified the latent random-effects Hessian as tridiagonal:

- Random effects: 20
- Hessian nonzeros: 58
- Detected structure: tridiagonal
- Solver complexity: O(n)

This enabled use of the specialized tridiagonal Laplace backend instead of a general sparse or dense factorization.

### Performance

| Metric | Quadra | TMB |
|----------|----------|----------|
| End-to-end runtime | ~0.18 s | ~0.49 s |
| Peak RSS | ~5 MB | ~170 MB |

In this benchmark Quadra achieved equivalent results while using substantially less memory and lower end-to-end runtime through automatic structure detection and specialized Laplace evaluation.
