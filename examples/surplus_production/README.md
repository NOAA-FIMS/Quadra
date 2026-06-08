# Surplus production example

This is a deterministic Schaefer surplus production model:

```text
B[t+1] = B[t] + r B[t] (1 - B[t] / K) - C[t]
I[t]   = q B[t] exp(epsilon[t])
```

The example reports the objective, biomass trajectory, fitted index, residuals, and reference points:

```text
MSY = rK / 4
B_MSY = K / 2
F_MSY = r / 2
terminal depletion = B_terminal / K
```

Run from the repository root:

```bash
./run_surplus_production_example.sh
```

Suggested next steps:

1. Add finite-difference gradient checks.
2. Add a simple optimizer wrapper.
3. Add process error.
4. Add random effects and Laplace evaluation.

## Fit the demo data

A simple finite-difference optimizer example is available:

```bash
./run_fit_surplus_production_example.sh
```

This estimates:

```text
log_r
log_K
log_q
log_sigma_index
logit_B0_frac
```

The optimizer is intentionally simple. It is meant to establish the fisheries-facing
parameter-estimation workflow before replacing finite-difference gradients with
Quadra automatic differentiation.
