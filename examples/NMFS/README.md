# NMFS Assessment Examples

This directory contains synthetic fisheries stock-assessment examples
implemented with Quadra. They demonstrate assessment-shaped workflows; none is
an official NMFS assessment or management product.

These examples are application-oriented and are separated from smaller framework
examples so that the repository clearly distinguishes between:

- core Quadra demonstrations
- fisheries assessment model applications
- validation and comparison studies

## Current examples

| Example | Primary command | Scope |
| --- | --- | --- |
| AFSC Walleye Pollock | `./examples/NMFS/afsc_walleye_pollock/run_walleye_pollock_example.sh` | Compact assessment and scaling diagnostics |
| SEFSC Red Snapper | `./examples/NMFS/sefsc_red_snapper/run_red_snapper_quadra_fit.sh` | Laplace fit and TMB comparison |
| PIFSC Opakapaka | `./examples/NMFS/pifsc_opakapaka/run_opakapaka_projection.sh` | Fit, uncertainty, and projections |
| PIFSC Spatial Tuna | `make -C examples/NMFS/pifsc_spatial_tuna test-fast` | End-to-end spatial assessment |

## Prerequisites and conventions

The Quadra-only examples require Bash, a C++17 compiler, and the vendored
dependencies in `external/`. Run the commands from a complete repository
checkout. The shell runners locate the repository root automatically, so they
also work when invoked from their own example directory.

Set `CXX` to select a compiler, for example:

```bash
CXX=g++ ./examples/NMFS/pifsc_opakapaka/run_opakapaka_projection.sh
```

The optional TMB comparisons require R, the R package `TMB`, and a working R
C++ toolchain. Spatial Tuna reporting additionally requires `Rscript` and the R
packages documented in its README. Generated executables are placed under
`build/`; assessment products are written below each example's `outputs/` or
local `build/` directory.

The default commands below use synthetic inputs included in the repository.
Do not substitute confidential or controlled assessment data into a public
checkout.

## Quick smoke check

These commands cover the short public entry points:

```bash
./examples/NMFS/afsc_walleye_pollock/run_walleye_pollock_example.sh
./examples/NMFS/sefsc_red_snapper/run_red_snapper_level0.sh
./examples/NMFS/pifsc_opakapaka/run_opakapaka_projection.sh
make -C examples/NMFS/pifsc_spatial_tuna test-fast
```

The full fits and simulation-recovery workflows can take substantially longer;
use the per-example README before launching them.

### AFSC Walleye Pollock

Path: `examples/NMFS/afsc_walleye_pollock`

This example includes catch, index, and age-composition observations plus an
optional recruitment-deviation scaling diagnostic.

### SEFSC Red Snapper

Path: `examples/NMFS/sefsc_red_snapper`

This example includes:

- age-structured population dynamics
- recruitment deviations as random effects
- Laplace approximation
- exact gradient validation
- comparison against a TMB implementation

The Red Snapper example is currently treated as a completed validation model for
Quadra's exact Laplace machinery.

### PIFSC Opakapaka

Path: `examples/NMFS/pifsc_opakapaka`

This example includes:

- Pacific Islands assessment-style projection workflow
- synthetic data input
- uncertainty reporting
- derived quantities
- projection uncertainty outputs
- comparison against a TMB implementation

### Spatial Tuna Assessment

Path: `examples/NMFS/pifsc_spatial_tuna`

This synthetic, end-to-end assessment example includes:

- seasonal age-structured, multi-region, multi-fleet population dynamics
- recruitment deviations marginalized with Quadra's structure-aware Laplace machinery
- phased fitting, simulation recovery, diagnostics, and runtime configuration
- exact posterior sampling with pCN, AD-NUTS, and native nonlinear transport proposals
- conditional random-effect reconstruction, reference points, and projections
- a detailed mathematical and computational specification in `MODEL_AND_QUADRA.md`
