# State-space surplus production example

This example extends the deterministic Schaefer surplus production model with
log-scale process deviations in biomass dynamics.

The joint objective is:

```text
joint_nll(theta, u)
```

where:

```text
theta = fixed effects
u     = annual log-scale process deviations
```

The process model is:

```text
pred_B[t+1] = B[t] + r B[t] (1 - B[t] / K) - C[t]
log_B[t+1]  = log(pred_B[t+1]) + u[t]
u[t]        ~ Normal(0, sigma_process)
```

The observation model is:

```text
log(I[t]) ~ Normal(log(q) + log_B[t], sigma_index)
```

Run:

```bash
bash examples/state_space_surplus_production/run_state_space_surplus_production_joint_example.sh
```

Next steps:

1. Add finite-difference gradient checks for `joint_objective(theta, u)`.
2. Add a Newton solver for `u_hat(theta)`.
3. Add the Laplace correction.
4. Replace finite differences with Quadra exact gradients.

## Random-effects-only fit

Hold fixed effects constant and optimize the latent process deviations:

```bash
bash examples/state_space_surplus_production/run_fit_state_space_surplus_u_example.sh
```

This computes:

```text
u_hat(theta) = argmin_u joint_nll(theta, u)
```

This is the step immediately before adding the Laplace correction.

## Dense Laplace example

Evaluate the Laplace approximation for the state-space surplus production model:

```bash
bash examples/state_space_surplus_production/run_state_space_surplus_dense_laplace_example.sh
```

This version uses:

```text
finite-difference gradient for u optimization
finite-difference dense H_uu
dense Cholesky log determinant
```

It is intentionally transparent. The next step is to replace the finite-difference
pieces with Quadra exact derivatives and then compare the same model to TMB.
