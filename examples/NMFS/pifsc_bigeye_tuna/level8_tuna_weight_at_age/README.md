# Level 8: Tuna Weight-at-Age

Level 6 had the best structural support so far, but recruitment diagnostics
still showed persistent recruitment deviations. Inspection of the shared
age-structured header showed that bigeye tuna was being modeled with a very
small red-snapper-like weight-at-age schedule.

This level keeps the Level 6 model structure and changes only weight-at-age to
a tuna-like schedule.

## Model change

Compared with Level 6:

- keep purse-seine fixed age-based selectivity
- keep longline logistic selectivity
- keep fleet-specific q
- keep annual recruitment deviations
- replace toy/red-snapper-like weight-at-age with tuna-like weight-at-age

## Scientific hypothesis

If persistent recruitment deviations were compensating for unrealistically low
adult weight-at-age, then replacing the weight schedule should change biomass
scale, q, Fbar, and recruitment diagnostics.

## Decision metrics

Compare Level 6 vs Level 8:

- objective
- R0
- Fbar
- q values
- rec_prior_nll
- recruitment sd
- recruitment lag-1 correlation
- age-comp NLL
