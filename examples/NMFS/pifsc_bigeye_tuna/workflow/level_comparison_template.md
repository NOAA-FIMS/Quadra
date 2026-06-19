# Model Level Comparison Template

Use this template every time a new level is added.

## Level name

Example:

```text
Level 1B: explicit fleet catchability
```

## Model change

Describe the single feature added.

## Scientific hypothesis

What should this feature explain or identify?

## Data expected to identify the feature

List the data source(s):

- catch
- index
- composition
- tagging
- priors

## Fit diagnostics

```text
objective:
gradient_norm:
converged:
max_gradient_parameter:
max_gradient_value:
iterations:
message:
```

## Random-effect Hessian diagnostics

```text
positive_definite:
condition_number:
min_eigenvalue:
max_eigenvalue:
structural_nonzeros:
effective_bandwidth_90:
effective_sparsity_90:
max_abs_correlation:
```

## Functional analysis interpretation

Top influential parameters:

```text
1.
2.
3.
```

Top correlation pairs:

```text
1.
2.
3.
```

## Identifiability interpretation

Weak or suspicious directions:

```text
-
-
-
```

## Scientific interpretation

What did the diagnostics teach us?

## Decision

Choose one:

- retain feature
- simplify feature
- pool feature
- regularize feature
- remove feature
- collect/add more data before estimating feature

## Next step

What is the next single model change?
