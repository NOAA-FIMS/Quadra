# Level 13: Estimated Initial Numbers

Level 12 corrected inherited reef-fish life history and converged cleanly, but
the recruitment deviations still showed a persistent positive block in the early
years. That suggests recruitment may be repairing a bad initial age structure.

Level 13 keeps the Level 12 assumptions:

- fixed M = 0.18
- fixed longline q = 0.00005
- tuna-like provisional weight-at-age
- tuna-like provisional maturity-at-age

and adds fixed-effect initial log-number deviations by age:

```text
init_log_number_dev_age_1 ... init_log_number_dev_age_10
```

These deviations multiply the equilibrium initial numbers-at-age before the
first modeled year.

A normal prior is applied to each initial-number deviation so this is a
diagnostic flexibility experiment, not an unconstrained free initial state.

## Scientific question

Does estimating initial numbers reduce the persistent early recruitment block?
