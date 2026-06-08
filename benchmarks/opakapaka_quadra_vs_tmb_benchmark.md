# Quadra vs TMB: Synthetic Opakapaka-Style Projection Benchmark

Synthetic and public-data-safe. Not an official assessment.

## Summary

| Engine | Objective | Wall time ms | Peak RSS MB | Grad norm | Converged | Random effects | Structure | Backend |
|---|---:|---:|---:|---:|---|---:|---|---|
| Quadra | -58.679607 | 140.000 | 14.922 | 0.000045 | yes | 30 | tridiagonal | tridiagonal |
| TMB | -26.971822 | 550.000 | 190.656 | — | yes | 30 | — | — |

## Observed Ratios

- Quadra was approximately **3.9× faster** end-to-end.
- Quadra used approximately **12.8× less peak memory**.
- Quadra detected the random-effect Hessian as **tridiagonal** with **88 nonzeros**.
- Quadra used a specialized **O(n) tridiagonal Laplace backend**.

## Convergence

The Quadra optimizer reached the requested fixed-effect gradient tolerance at outer evaluation 9:

```text
L-BFGS: outer eval = 9
fx        = -58.679607
|grad|    = 0.000045
```

A follow-up optimizer patch now stops at the first iterate satisfying the requested gradient tolerance rather than continuing into later line-search instability.

## Significance

This benchmark demonstrates that Quadra's structure-aware Laplace implementation can deliver substantial runtime and memory reductions for mixed-effects projection workflows while automatically detecting exploitable Hessian structure.

## Next Benchmark

Scale the random-effect dimension while keeping the model formulation fixed:

- 30 random effects
- 60 random effects
- 120 random effects
- 240 random effects
- 480 random effects
- 960 random effects

Record:

- End-to-end runtime
- Peak RSS memory
- Hessian nonzeros
- Structure classification
- Backend selected

This experiment will quantify the scaling behavior of Quadra's structure-aware Laplace engine relative to TMB.
