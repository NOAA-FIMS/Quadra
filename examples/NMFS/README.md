# NMFS Assessment Examples

This directory contains fisheries stock assessment examples implemented with
Quadra.

These examples are application-oriented and are separated from smaller framework
examples so that the repository clearly distinguishes between:

- core Quadra demonstrations
- fisheries assessment model applications
- validation and comparison studies

## Current examples

### SEFSC Red Snapper

Path:

```text
examples/NMFS/sefsc_red_snapper
```

This example includes:

- age-structured population dynamics
- recruitment deviations as random effects
- Laplace approximation
- exact gradient validation
- comparison against a TMB implementation

The Red Snapper example is currently treated as a completed validation model for
Quadra's exact Laplace machinery.

### PIFSC Opakapaka

Path:

```text
examples/NMFS/pifsc_opakapaka
```

This example includes:

- Pacific Islands assessment-style projection workflow
- synthetic data input
- uncertainty reporting
- derived quantities
- projection uncertainty outputs
- comparison against a TMB implementation

### Spatial Tuna Assessment

Path:

```text
examples/NMFS/pifsc_spatial_tuna
```

This synthetic, end-to-end assessment example includes:

- seasonal age-structured, multi-region, multi-fleet population dynamics
- recruitment deviations marginalized with Quadra's structure-aware Laplace machinery
- phased fitting, simulation recovery, diagnostics, and runtime configuration
- exact posterior sampling with pCN, AD-NUTS, and native nonlinear transport proposals
- conditional random-effect reconstruction, reference points, and projections
- a detailed mathematical and computational specification in `MODEL_AND_QUADRA.md`
