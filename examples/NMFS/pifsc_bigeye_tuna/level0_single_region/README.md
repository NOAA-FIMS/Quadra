# Level 0: Single-Region Bigeye Tuna Prototype

This level is the baseline model for the Bigeye diagnostic-guided modeling
ladder.

It intentionally starts simple.

## Purpose

Create a tuna-like single-region model with the same diagnostic package used by
the Red Snapper example.

## Initial model

- one region
- one aggregate catch series
- one abundance index
- age-structured dynamics
- recruitment deviations as random effects
- fixed effects for recruitment scale, fishing mortality, and catchability

## Key question

Can the available data separate recruitment scale, fishing mortality,
catchability, and recruitment deviations?

## Expected outputs

- `bigeye_fit_summary.csv`
- `bigeye_objective_components.csv`
- `bigeye_functional_analysis_report.txt`
- `bigeye_functional_analysis_report.csv`
- `bigeye_laplace_structure_report.txt`
- `bigeye_laplace_structure_report.csv`
- `bigeye_reference_points.csv`
