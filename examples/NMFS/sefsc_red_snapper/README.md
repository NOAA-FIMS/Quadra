# SEFSC Red-Snapper-Style Assessment Example

This is a synthetic, public-data-safe red-snapper-style assessment example. It
is not an official assessment and should not be used for management advice.

It provides an SEFSC-style validation case for Quadra with age structure,
selectivity, recruitment deviations, diagnostics, and a side-by-side TMB
implementation.

## Implemented model features

- age-structured population dynamics
- catch likelihood
- survey/index likelihood
- age-composition likelihood
- recruitment deviations as random effects
- age-based selectivity
- fitted trajectories and residual diagnostics
- objective-component, Hessian-structure, and functional-analysis reports
- reference-point summaries

## Requirements

Quadra-only runs require Bash and a C++17 compiler. The repository vendors
Eigen and LBFGSpp. Set `CXX` to select a compiler. The optional comparison also
requires `Rscript`, the R package `TMB`, Python 3, and a working R C++ toolchain.

## Build and run

The runners work from anywhere inside the repository checkout. Start with the
fast deterministic stages:

```bash
./examples/NMFS/sefsc_red_snapper/run_red_snapper_level0.sh
./examples/NMFS/sefsc_red_snapper/run_red_snapper_objective.sh
./examples/NMFS/sefsc_red_snapper/run_red_snapper_age_structured.sh
```

Run the full Quadra Laplace fit with:

```bash
./examples/NMFS/sefsc_red_snapper/run_red_snapper_quadra_fit.sh
```

Run the optional Quadra/TMB comparison with:

```bash
./examples/NMFS/sefsc_red_snapper/run_quadra_vs_tmb_comparison.sh
```

Executables are written to `build/examples/`. CSV and text reports are written
to `examples/NMFS/sefsc_red_snapper/outputs/`; repeated runs replace files with
the same names. Successful fit output reports `converged: yes` and exits zero.

## Directory layout

```text
data/        synthetic or public-data-safe inputs
quadra/      Quadra implementation
tmb/         TMB reference implementation
outputs/     generated and reference outputs
validation/  comparison summaries and validation notes
```

Some outputs are versioned as reference artifacts. Review diffs before
committing a run made with a different compiler or configuration.

## Validation scope

The current validation workflow covers:

1. deterministic age-structured dynamics,
2. one abundance index,
3. synthetic catch observations,
4. recruitment deviations,
5. TMB side-by-side comparison,
6. structure and functional diagnostics.

The TMB comparison checks numerical agreement for selected fit fields; it is
not a claim that either synthetic model reproduces an operational assessment.
