# Diagnostic-Guided Model Construction

The goal is to make model development auditable.

For each model level, record:

## Model change

What was added?

Examples:

- new fleet
- fleet-specific selectivity
- length composition
- regional state structure
- movement
- tagging likelihood

## Expected information gain

What should the new feature identify?

Examples:

- catchability
- fleet selectivity
- recruitment variability
- regional abundance
- movement probability

## Diagnostics to compare

- convergence status
- fixed-effect gradient norm
- Huu positive definiteness
- Huu condition number
- effective sparsity
- effective bandwidth
- random-effect correlation structure
- top parameter influence
- weakest curvature directions
- objective component contributions

## Decision

Retain, simplify, or remove the feature.

A feature should not be retained only because it is biologically plausible. It
should also improve the model in a way that the data and diagnostics support.
