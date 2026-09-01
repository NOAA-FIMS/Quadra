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

Requirements are Bash, a C++17 compiler, and the vendored dependencies in the
repository. From the repository root:

```bash
./examples/NMFS/pifsc_opakapaka/run_opakapaka_projection.sh
```

The same script can be invoked as `./run_opakapaka_projection.sh` from this
directory. Set `CXX` to choose a compiler. The executable is written to
`build/examples/pifsc_opakapaka`; generated CSV files go to this example's
`outputs/` directory.

Outputs are written to:

```text
examples/NMFS/pifsc_opakapaka/outputs/synthetic_fit_summary.csv
examples/NMFS/pifsc_opakapaka/outputs/synthetic_projection_scenarios.csv
```

## Opakapaka projection validation

A synthetic opakapaka-style state-space projection model was implemented in
both Quadra and TMB using identical data and model structure. The values below
are a recorded validation baseline, not assertions about every compiler,
machine, or later optimizer configuration. Current-run values in `outputs/`
are authoritative for the checked-out revision.

### Numerical Agreement

| Metric | Quadra | TMB |
|----------|----------|----------|
| Objective | 38.888259 | 38.888010 |
| log(q) | -6.968967 | -6.964753 |
| q | 0.000941 | 0.000945 |
| r | 0.388383 | 0.386970 |
| K | 1337.136918 | 1334.502402 |
| Random effects | 20 | 20 |

The two implementations converge to closely matching parameter estimates and
objective values while using their respective Laplace-gradient implementations.

### Quadra Structure Detection

Quadra automatically identified the latent random-effects Hessian as tridiagonal:

- Random effects: 20
- Hessian nonzeros: 58
- Detected structure: tridiagonal
- Solver complexity: O(n)

This enabled use of the specialized tridiagonal Laplace backend instead of a general sparse or dense factorization.

### Recorded performance

| Metric | Quadra | TMB |
|----------|----------|----------|
| End-to-end runtime | ~0.18 s | ~0.49 s |
| Peak RSS | ~5 MB | ~170 MB |

These timings and memory measurements are illustrative and are not portable
benchmarks. Compiler, optimization flags, hardware, and warm/cold process state
all affect them. In the recorded run, Quadra achieved equivalent numerical
results while using structure detection and a specialized Laplace evaluation.
