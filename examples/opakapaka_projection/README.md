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
examples/opakapaka_projection/outputs/synthetic_fit_summary.csv
examples/opakapaka_projection/outputs/synthetic_projection_scenarios.csv
```
