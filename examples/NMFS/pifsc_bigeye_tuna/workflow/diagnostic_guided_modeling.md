# Diagnostic-Guided Model Construction

Quadra's purpose in this prototype is not only to fit a model. The purpose is to
make model development auditable.

## Workflow

1. Start with the simplest defensible model.
2. Fit the model.
3. Generate diagnostics.
4. Add one model feature.
5. Refit.
6. Compare diagnostics.
7. Retain, simplify, or remove the feature.

## Required outputs at each level

- fit summary
- objective components
- functional analysis report
- Laplace structure report
- reference points, when biologically meaningful
- model decision notes

## Future diagnostic frontier: identifiability report

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
