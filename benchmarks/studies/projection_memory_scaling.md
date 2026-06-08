# Projection Memory Scaling Study

This benchmark study evaluates memory and runtime behavior for long-horizon latent-state projection workloads.

The goal is to measure whether a structure-aware Laplace runtime can keep memory requirements modest as the number of latent states increases. This is especially relevant for assessment workflows that combine long historical time series with long projection horizons.

The benchmark is intentionally framed in terms of computational behavior rather than a specific application interface. The comparison uses:

- **Structure-aware runtime:** Quadra `LaplaceEvaluator` using persistent latent-state storage, warm-start Newton updates, direct tridiagonal Hessian values, and a specialized tridiagonal log-determinant path.
- **Generic AD/Laplace runtime:** the current TMB benchmark implementation used for this state-space surplus production test.

## Benchmark model

The benchmark uses a state-space surplus production model with latent biomass states. The random-effects Hessian for this model is tridiagonal, corresponding to first-order Markov dependence among latent states.

A model-analysis report for this Hessian identifies:

- detected structure: tridiagonal
- recommended solver: Newton
- recommended backend: tridiagonal
- expected complexity: `O(n)`

This structure is what allows the optimized runtime to avoid generic sparse Hessian rebuilds and generic sparse factorizations.

## Projection scenarios

The scenarios below are expressed as historical years plus projection years. The latent-state count is represented as the sum of those two quantities for this benchmark.

| History years | Projection years | Latent states | Structure-aware RSS | Structure-aware runtime | Generic AD/Laplace RSS | Generic AD/Laplace runtime | RSS reduction |
|--------------:|-----------------:|--------------:|--------------------:|------------------------:|-----------------------:|---------------------------:|--------------:|
| 30 | 30 | 60 | 1.3 MB | 0.006 ms | 191.5 MB | 0.800 ms | 148× |
| 100 | 100 | 200 | 1.3 MB | 0.025 ms | 313.5 MB | 9.100 ms | 242× |
| 300 | 100 | 400 | 1.4 MB | 0.034 ms | 663.2 MB | 40.400 ms | 488× |
| 300 | 300 | 600 | 1.5 MB | 0.038 ms | 1.29 GB | 101.900 ms | 862× |
| 1000 | 1000 | 2000 | 1.9 MB | 0.135 ms | 8.07 GB | 1947.000 ms | 4,266× |
| 10000 | 10000 | 20000 | 6.6 MB | 1.291 ms | — | — | — |
| 50000 | 50000 | 100000 | 57.6 MB | 6.645 ms | — | — | — |

## Key observations

At 2,000 latent states, the structure-aware runtime used 1.9 MB peak RSS, while the generic AD/Laplace implementation used 8.07 GB. That is a measured RSS reduction of approximately 4,266×.

The largest stress case in this study used 100,000 latent states and completed with 57.6 MB peak RSS and 6.645 ms warm-start runtime.

Across the measured projection scenarios, the structure-aware runtime remained well below a 16 GB continuous-integration memory limit. The largest measured stress case, 100,000 latent states, used less than 100 MB peak RSS in this benchmark run.

The generic AD/Laplace implementation showed substantially higher memory use in the small and moderate scenarios that were run back-to-back. Larger generic-runtime scenarios were skipped to avoid long runtimes and excessive memory use.

## Why this matters

Long projection horizons can create large latent-state vectors. If the computational backend treats the model as a generic sparse problem, memory use can grow quickly enough to become a limiting factor for testing, continuous integration, and routine assessment workflows.

The results here suggest that when latent-state structure is known or detectable, a structure-aware runtime can substantially reduce peak memory while preserving very small per-evaluation runtimes.

## Interpretation

These results should be interpreted as a benchmark of one model class: a tridiagonal state-space surplus production model. They do not imply that all mixed-effects models will show the same memory reduction.

The strongest conclusion is narrower and more defensible:

> For projection-heavy latent-state models with exploitable tridiagonal structure, a structure-aware Laplace runtime can reduce memory requirements by orders of magnitude relative to a generic AD/Laplace workflow.

## Reproducing the study

The benchmark was generated with:

```bash
TMB_MAX_N=2000 ./run_fims_projection_memory_suite.sh 10
```

The summary table used for this report was:

```text
benchmarks/fims_projection_memory_suite_20260604_142643/summary.tsv
```

For larger structure-aware stress tests:

```bash
./run_state_space_laplace_memory_benchmark.sh 10 10000,20000,50000,100000
```

## Limitations and next benchmarks

Recommended follow-up benchmarks:

1. Measure peak RSS during full outer optimization, not only fixed-theta evaluation.
2. Add AR(1) and random-walk recruitment examples.
3. Add random-intercept and block-diagonal examples.
4. Add a theta-perturbation benchmark that reuses warm starts during outer optimization.
5. Compare generic-runtime memory under carefully bounded large-state scenarios.
