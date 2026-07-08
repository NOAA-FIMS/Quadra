# Level 7: Dual Age-Based Selectivity

Level 6 replaced purse-seine logistic selectivity with fixed age-based
selectivity. That substantially improved the fit and returned R0/q/Fbar to
more sensible scales. Recruitment diagnostics then showed a persistent,
smooth recruitment trajectory, suggesting that recruitment may still be
absorbing remaining structural mismatch.

Level 7 tests whether the remaining mismatch is caused by longline logistic
selectivity.

## Model change

Compared with Level 6:

- replace longline logistic selectivity with fixed age-based selectivity
- keep purse-seine fixed age-based selectivity
- estimate only:
  - log_r0
  - log_fbar
  - log_q_longline
  - log_q_purse_seine
- retain annual recruitment deviations as random effects
- retain Level 6 recruitment diagnostics

## Scientific hypothesis

If Level 6 recruitment deviations were compensating for remaining longline
selectivity misspecification, then Level 7 should reduce recruitment
persistence and recruitment prior burden.

## Decision metrics

Compare Level 6 vs Level 7:

- objective
- age_comp_nll
- rec_prior_nll
- recruitment sd
- recruitment lag-1 correlation
- recruitment roughness
- fixed-effect gradient
- safe wiggle behavior

## Interpretation

If recruitment persistence drops, then recruitment was compensating for
selectivity structure. If recruitment persistence remains high, the persistent
recruitment signal may reflect data construction, catch/index weighting, or
a need for an explicit recruitment process rather than independent annual
random effects.
