# State-space surplus production scaling benchmark

This benchmark compares fixed-theta Laplace objective evaluation and isolated
process peak resident memory (RSS) for the same
state-space surplus production model in:

```text
Quadra persistent latent-state tridiagonal implementation
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

| n | Quadra ms | TMB ms | Speedup | Quadra RSS MiB | TMB RSS MiB | RSS ratio |
|---:|---:|---:|---:|---:|---:|---:|
| 25 | 0.011 | 0.2 | 18.2x | 1.09 | 199.52 | 182.4x |
| 50 | 0.026 | 0.7 | 27.4x | 1.06 | 192.00 | 180.7x |
| 100 | 0.052 | 2.4 | 46.2x | 1.11 | 217.28 | 195.9x |
| 250 | 0.113 | 16.0 | 141.1x | 1.16 | 398.25 | 344.4x |
| 500 | 0.228 | 73.9 | 324.1x | 1.39 | 1000.58 | 719.5x |
| 1000 | 0.423 | 390.7 | 923.8x | 1.88 | 2717.70 | 1449.4x |

Objectives matched to numerical precision.

## Interpretation

This benchmark does **not** show that Quadra is universally faster than TMB.

It shows that when the model has known latent Markov structure and Quadra
exploits persistent tridiagonal structure, the Laplace evaluation and its
memory footprint can scale far better than a generic AD/Laplace path.

## Run

From the repository root:

```bash
./benchmarks/state_space_surplus_scaling/run_benchmark.sh \
  10 25,50,100,250,500,1000
```

Outputs are written to:

```text
benchmarks/state_space_surplus_scaling/results.csv
benchmarks/state_space_surplus_scaling/scaling_plot.png
```

Each model and size runs in a separate measured process. Compilation happens
before measurement, and peak RSS is normalized to MiB in `results.csv`.
