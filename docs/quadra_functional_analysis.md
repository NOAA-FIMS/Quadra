# Quadra Functional Analysis v1

## Purpose

Quadra Functional Analysis v1 summarizes whether a mixed-effects model fit is
numerically healthy, how its latent-state curvature is structured, and whether
the apparent Hessian complexity is truly global or effectively local.

## v1 Diagnostics

- Model health assessment
- Optimization diagnostics
- Curvature diagnostics
- Spectral structure
- Huu structure
- Effective sparsity
- Effective bandwidth
- Uncertainty summary
- Parameter influence
- Parameter geometry
- Correlation graph topology
- Latent state summary
- Markdown and CSV reporting

## Pollock Showcase

The synthetic AFSC walleye-pollock-style example currently generates:

- `examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_analysis.md`
- `examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_functional_analysis_report.txt`
- `examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_functional_analysis_report.csv`
- `examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_laplace_structure_report.txt`
- `examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_laplace_structure_report.csv`

### Key Pollock Message

The Pollock example demonstrates why functional analysis is more informative
than symbolic sparsity alone: the random-effect Hessian can appear structurally
dense while the numerical curvature and uncertainty graph reveal local,
AR(1)-style dependence.

## Modernization Targets

### Opakapaka

Scaffolded at:

`examples/NMFS/pifsc_opakapaka`

Next target: wire the existing Quadra driver to the Functional Analysis v1
report API and generate a matching markdown report.

### Red Snapper

Scaffolded at:

`examples/NMFS/sefsc_red_snapper`

Next target: wire the existing Quadra driver to the Functional Analysis v1
report API and generate a matching markdown report.

## Meeting Talking Point

Functional Analysis v1 turns raw optimizer output into a reusable model-review
artifact. Instead of only reporting objective values and convergence codes, it
summarizes optimization health, curvature health, latent-state structure,
uncertainty topology, and numerical compressibility.
