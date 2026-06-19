# Diagnostic-Guided Model Construction

Quadra's purpose in the Bigeye prototype is not only to fit a model. The
purpose is to make model development auditable.

## Core workflow

1. Start with the simplest defensible model.
2. Fit the model.
3. Generate diagnostics.
4. Interpret diagnostics scientifically.
5. Add one model feature.
6. Refit.
7. Compare diagnostics.
8. Retain, simplify, pool, regularize, or remove the feature.

## Required outputs at each level

Each level should write:

- fit summary
- objective components
- functional analysis report
- Laplace structure report
- reference points, when biologically meaningful
- model decision notes
- scientific reasoning log entry

## Diagnostic categories

### Optimization health

- convergence status
- objective value
- gradient norm
- iterations
- optimizer message
- maximum fixed-effect gradient parameter

### Curvature health

- Huu positive definiteness
- minimum eigenvalue
- condition number
- effective rank
- high-correlation random-effect pairs

### Structural health

- structural nonzeros
- effective sparsity
- effective bandwidth
- compression relative to structural matrix

### Parameter influence

- variance share
- correlation centrality
- curvature column norm
- influence ranking

### Identifiability health

The next diagnostic frontier is identifying weak parameter combinations.

The goal is to expose weak directions such as:

```text
+ log_q_fleet_1
- abundance_region_1
```

or:

```text
+ movement_region_1_to_2
- recruitment_region_2
```

These directions indicate parameters that the data cannot cleanly separate.

## Scientific rule

A biologically plausible feature should not automatically be retained.

It should be retained only if the diagnostics show that the available data
support estimating it.
