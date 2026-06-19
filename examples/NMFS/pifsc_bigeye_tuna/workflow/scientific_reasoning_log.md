# Bigeye Tuna Diagnostic-Guided Scientific Reasoning Log

This document records the scientific reasoning behind each model-building step.

The goal is not to build a large model quickly. The goal is to build a model
whose complexity is justified by diagnostics.

## Working philosophy

A model feature should earn its place.

For each new feature, we ask:

1. What biological or observation-process hypothesis motivates the feature?
2. What data should identify the feature?
3. What diagnostics changed after adding the feature?
4. Did the feature improve interpretation, stability, or information content?
5. Did the feature introduce confounding, weak curvature, or overparameterization?
6. Should the feature be retained, simplified, pooled, regularized, or removed?

## Diagnostic interpretation rules

### Random-effect Hessian health

If the random-effect Hessian is positive definite and well-conditioned, but the
fixed-effect optimizer stalls, the first suspect is not the Laplace machinery.
The first suspect is fixed-effect scaling, data conflict, or weak
identifiability.

### Fixed-effect gradient health

If the maximum fixed-effect gradient is concentrated on one parameter, that
parameter becomes the first diagnostic target.

Example:

```text
max_gradient_parameter: log_q
```

This suggests catchability/abundance scaling may be weakly identified or
conflicted with other fixed effects.

### Scale confounding

If recruitment scale increases while catchability decreases, the model may be
trading abundance scale against observation scale.

Example pattern:

```text
R0 increases
q decreases
index fit remains similar
```

Interpretation:

The data may not separate absolute abundance from catchability.

### Complexity rule

Do not add movement, regions, tagging, or fleet-specific selectivity when the
current simpler model already shows unresolved q/F/R0 confounding.

## Level 0: single-region baseline

### Model change

The Level 0 model starts with:

- one region
- one aggregate catch series
- one abundance index
- age-structured dynamics
- recruitment deviations as random effects
- fixed effects for recruitment scale, fishing mortality, catchability, and
  selectivity

### Scientific purpose

Establish a baseline diagnostic package before adding fleet structure.

### Diagnostic expectation

Level 0 should answer whether the base model is numerically healthy and whether
the random effects are identifiable under a simple observation process.

### Decision

Level 0 is the baseline. It should not be treated as a tuna assessment. It is a
diagnostic control case.

## Level 1A: aggregated multi-fleet scaffold

### Model change

Two synthetic fleets were introduced in the data file:

- longline
- purse seine

The first scaffold aggregated fleet catches by year and averaged fleet indices
by year so the existing Level 0 biological model could run unchanged.

### Scientific hypothesis

If fleet information is aggregated into one catch series and one index, the
model should behave similarly to Level 0 unless the aggregation creates a scale
conflict.

### Diagnostics observed

The first Level 1A scaffold produced:

```text
objective:                 115.031861
gradient_norm:             0.429520
converged:                 no
max_gradient_parameter:    log_q
max_gradient_value:        0.346158
Huu positive definite:     yes
Huu condition number:      6.20
structural nonzeros:       362 / 400
max random-effect corr:    about 0.21
```

### Scientific interpretation

The random-effect curvature remained healthy. The random-effect Hessian was
positive definite and well-conditioned.

The fixed-effect optimizer stalled, with the largest remaining gradient on
`log_q`.

The aggregated data created a likely scale conflict:

```text
catch = longline catch + purse seine catch
index = average(longline index, purse seine index)
```

This roughly doubled removals while leaving the index scale comparable to a
single-fleet index.

The model responded by increasing the recruitment/abundance scale and reducing
catchability:

```text
R0 increased
q decreased
```

This is consistent with abundance/catchability scale confounding.

### Decision

Level 1A is a useful diagnostic failure.

It should not be treated as a validated multi-fleet model.

The next step is not to add regions, movement, or tagging. The next step is to
make the fleet observation process explicit.

## Level 1B: explicit fleet catchability

### Planned model change

Introduce explicit fleet-level observation terms:

- shared biological population
- fleet-specific catch/index observations
- fleet-specific catchability parameters, initially:
  - `log_q_longline`
  - `log_q_purse_seine`

### Scientific hypothesis

If Level 1A failed because of observation-process aggregation, then allowing
fleet-specific catchability should reduce the `log_q` scaling conflict.

### Diagnostic questions

- Does the model converge?
- Does the maximum fixed-effect gradient decrease?
- Does the maximum gradient remain on a catchability parameter?
- Do fleet-specific q values become identifiable?
- Does Huu remain positive definite?
- Does fixed-effect correlation reveal q/R0/F confounding?
- Does the functional analysis show different influence patterns than Level 1A?

### Decision rule

If explicit fleet q improves convergence and reduces the log_q gradient, then
the problem was likely aggregation-induced observation-process confounding.

If explicit fleet q does not improve convergence, then the issue is more
structural and we should inspect fixed-effect Hessian geometry before adding
additional biological complexity.

## Level 2: fleet-specific selectivity

### Planned model change

Add fleet-specific selectivity only after Level 1B diagnostics show that
fleet-specific catchability is not enough.

### Scientific hypothesis

Different fleets may observe different parts of the age/size distribution.

### Diagnostic questions

- Does fleet-specific selectivity improve fit enough to justify extra
  parameters?
- Which fleet identifies selectivity?
- Are selectivity and q confounded?
- Are selectivity and F confounded?

### Decision rule

Fleet-specific selectivity should be simplified or pooled if diagnostics show
weak curvature or parameter redundancy.

## Level 3: spatial structure

### Planned model change

Add regional abundance states only after fleet observation processes are
diagnostically stable.

### Scientific hypothesis

Regional structure may explain differences among fleets and indices.

### Diagnostic questions

- Are regional abundance states separately identifiable?
- Does q become region-confounded?
- Does the Huu correlation graph show strong regional blocks?

### Decision rule

Do not add movement until regional abundance is at least partially identifiable.

## Level 4: movement

### Planned model change

Add movement among regions.

### Scientific hypothesis

Movement may explain spatial redistribution of abundance.

### Diagnostic questions

- Is movement identifiable separately from regional recruitment?
- Is movement identifiable separately from regional q?
- Does movement reduce residual structure or only absorb noise?

### Decision rule

Movement should be pooled, simplified, or regularized if weak directions show
movement/recruitment/q confounding.

## Level 5: tagging

### Planned model change

Add tagging data only after movement structure has a clear diagnostic need.

### Scientific hypothesis

Tagging may identify movement and reporting rates.

### Diagnostic questions

- Does tagging reduce movement uncertainty?
- Are movement and reporting rate separable?
- Which tag groups contribute information?

### Decision rule

Retain tagging complexity only if it materially improves movement
identifiability.
