# Level 6: Purse-Seine Age-Based Selectivity

Level 5 showed that strong fleet composition contrast substantially reduced
the broad R0/F/q confounding, but the remaining weakest direction was dominated
by `logit_sel_a50_purse_seine`.

That suggests the two-parameter logistic curve is the wrong structure for purse
seine. The purse-seine fleet mainly observes young fish, so a free logistic a50
near the lower age boundary is weakly identified.

## Model change

Compared with Level 5:

- keep longline logistic selectivity
- replace purse-seine logistic selectivity with fixed age-based selectivity
- remove purse-seine a50 and slope fixed effects
- keep fleet-specific q
- keep strong fleet composition contrast

## Scientific hypothesis

If the remaining weakness is caused by an inappropriate purse-seine logistic
parameterization, replacing it with age-based selectivity should:

- remove the weak purse-seine a50 direction
- increase the smallest fixed-effect eigenvalue
- reduce condition number
- preserve or improve optimizer behavior

## Decision rule

If geometry improves, the diagnostic-guided conclusion is that purse-seine
selectivity should be represented as an age-based observation process rather
than a freely estimated logistic curve in this synthetic case.
