# Level 1B: Fleet-Specific Catchability Bigeye Prototype

This level tests the first explicit observation-process fix after the Level 1A
aggregated multi-fleet diagnostic failure.

## Scientific hypothesis

Level 1A aggregated two fleet catches but averaged fleet indices. Diagnostics
showed a fixed-effect confounding triangle among:

- `log_r0`
- `log_fbar`
- `log_q`

Level 1B tests whether explicit fleet-specific catchability reduces this
confounding.

## Model change

Compared with Level 1A:

- keep one shared biological population
- keep one shared recruitment-deviation random-effect vector
- keep fleet-specific observations in the data
- replace one aggregate catchability parameter with:
  - `log_q_longline`
  - `log_q_purse_seine`

## Diagnostic questions

- Does convergence improve?
- Does maximum fixed-effect gradient decrease?
- Does the maximum gradient remain on q?
- Does fixed-effect geometry still show strong R0/F/q confounding?
- Are the two fleet q values separately identifiable?
- Does Huu remain positive definite?

## Decision rule

If fleet-specific q improves convergence and reduces the weak R0/F/q direction,
then Level 1A's failure was likely caused by aggregation-induced observation
scale conflict.

If not, the confounding is structural and should be addressed before adding
fleet-specific selectivity, regions, movement, or tagging.
