# State-space surplus production scaling benchmark

This benchmark compares fixed-theta Laplace objective evaluation for the same
state-space surplus production model in:

```text
Quadra analytic latent-state tridiagonal implementation
TMB AD/Laplace implementation
```

The model is:

```text
pred_B[t+1] = B[t] + r B[t] (1 - B[t] / K) - C[t]
log_B[t+1]  = log(pred_B[t+1]) + process error
log(I[t])   = log(q) + log_B[t] + observation error
```

The fixed effects are held constant:

```text
r              = 0.5
K              = 700
q              = 0.0024
sigma_process  = 0.15
sigma_index    = 0.10
B0/K           = 0.90
```

The benchmark evaluates the marginal negative log likelihood after integrating
over latent log-biomass states.

## Current result

Representative run:

| n | Quadra ms/eval | TMB ms/eval | Quadra speedup |
|---:|---:|---:|---:|
| 25 | 0.027 | 0.100 | 3.7x |
| 50 | 0.060 | 0.600 | 9.9x |
| 100 | 0.117 | 2.200 | 18.7x |
| 250 | 0.255 | 14.100 | 55.2x |

Objectives matched to numerical precision.

## Interpretation

This benchmark does **not** show that Quadra is universally faster than TMB.

It shows that when the model has known latent Markov structure and Quadra
exploits the analytic tridiagonal Hessian, the Laplace evaluation can scale
nearly linearly and substantially outperform a generic AD/Laplace path.

## Run

From the repository root:

```bash
./run_state_space_surplus_scaling_artifact.sh 10 25,50,100,250
```

Outputs are written to:

```text
benchmarks/state_space_surplus_scaling/results.csv
benchmarks/state_space_surplus_scaling/scaling_plot.png
```
