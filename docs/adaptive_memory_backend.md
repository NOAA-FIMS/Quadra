# Adaptive memory backend exploration

The first experimental seam chooses the directional Hdot batch size from a
graph-memory budget. Existing callers remain unchanged; adaptive behavior is
opt-in through `ExactGradientWorkspace::ResizeDirectionalBatchAdaptive`.

The planner accounts for current graph-owned reserved memory and estimates the
fixed incremental storage for each directional lane. It then selects the
largest batch that fits the budget, subject to configured minimum and maximum
sizes. A minimum batch of one provides a bounded-progress fallback when the
recorded base graph already consumes the budget.

```cpp
quadra::laplace::AdaptiveDirectionalBatchOptions options;
options.memory_budget_bytes = 512u * 1024u * 1024u;
options.maximum_batch_size = 32;

auto plan = workspace.ResizeDirectionalBatchAdaptive(theta.size(), options);
```

`TraceTermsSelectedInverseAdaptive` is the corresponding streaming executor.
It seeds each chunk using global direction indices, propagates the local batch,
contracts trace terms immediately, and advances without retaining an Hdot
matrix for every direction. It reports the chosen plan, number of batches, and
peak tracked graph-owned bytes.

Dynamic second-order edge nodes are model-dependent, so process RSS should
remain the final enforcement signal in a production backend.

An exploratory benchmark is available as:

```bash
make benchmarks/benchmark_adaptive_directional_batch
./benchmarks/benchmark_adaptive_directional_batch 1000 32 64
```

Its arguments are random-effect dimension, direction count, and memory budget
in MiB. Output is CSV so budget sweeps can be collected directly.

The public catch-at-age evaluator currently scales memory through persistent
Hdot worker tapes rather than the directional-batch workspace. Its focused
comparison benchmark is:

```bash
make benchmarks/benchmark_big_catch_at_age_hdot_workers
./benchmarks/benchmark_big_catch_at_age_hdot_workers 1 3
```

Run worker counts in separate processes so `ru_maxrss` gives a comparable peak
for each configuration.

## Dense catch-at-age path

Catch-at-age currently detects a dense random-effect Hessian and therefore uses
the polarized third-order backend rather than Hdot worker tapes. Its dense trace
is now contracted without materializing either Hdot or the full inverse:

1. solve for one inverse-Hessian column;
2. generate the lower-triangular polarized Hdot entries for that column;
3. accumulate diagonal entries once and off-diagonal entries twice;
4. discard the inverse column and continue.

This reduces additional dense trace-contraction storage from `O(n_random^2)`
to `O(n_random)`. The sparse factorization itself may still have quadratic fill
for a genuinely dense Hessian, and derivative runtime remains quadratic.

The catch-at-age benchmark accepts a final `0` or `1` argument to compare the
legacy materialized and streamed dense trace paths. On the 30-year model with
10 repetitions, streaming changed mean exact-evaluation runtime from 400.601 ms
to 402.542 ms (+0.48%), reduced measured peak RSS from 14.734 to 14.547 MiB,
and preserved the gradient checksum exactly.

The fitted Hessian is fully dense. Its residual after removing bandwidth 5 has
only 0.4185% of the full Frobenius norm, but remains numerical rank 30/30. This
makes approximate banding worth investigating, but does not support an exact
banded-plus-low-rank decomposition at strict tolerances.

## Banded approximation sweep

The fitted catch-at-age sweep applies the same bandwidth to the Hessian and its
directional derivative, making third-derivative work `O(n * bandwidth)` rather
than `O(n^2)`. With five repetitions, bandwidth 5 reduced gradient time from
413.918 ms to 141.151 ms (2.93x). Its log-determinant relative error was
5.7e-6 and mode-sensitivity error was 0.73%.

The gradient L2 error was 0.01093 while the reference gradient norm was only
0.00947. Thus bandwidth 5 is not accurate enough as a final convergence
gradient for this fitted case despite its strong objective and sensitivity
accuracy. It remains a candidate for early optimization, followed by the exact
streaming backend near convergence.

## Optimization trajectory replay

The exact catch-at-age trajectory contains 79 objective/gradient evaluations.
Bandwidth 5 was replayed at evaluations 1, 25, 50, 65, and 75. Gradient cosine
was at least 0.999999 through evaluation 50 (`|g| = 5.46`), then declined to
0.9715 at evaluation 65 (`|g| = 0.0435`) and 0.7072 near convergence.

A permanent switch to exact evaluation at the first `|g| <= 0.1` would assign
64 evaluations to bandwidth 5 and 15 to the exact backend. Using measured
per-gradient timings, this projects about a 52% reduction in derivative time.
The projection does not assume unchanged optimizer iterations; that must be
verified by the live hybrid-optimizer experiment. Exact objective line
searches and an exact final convergence check remain required safeguards.

## Live hybrid optimizer

The opt-in public hybrid optimizer uses one persistent evaluator, bandwidth 5
above a gradient norm of 0.1, and the full dense trace below it. The transition
is permanent. It recomputes the switch point exactly, preserves L-BFGS history
only while it continues to produce an exact descent direction, and requires an
exact convergence check.

On catch-at-age, exact-only converged in 67 iterations and 31.603 seconds. The
hybrid converged in 74 iterations using 54 approximate and 25 exact evaluations
and took 18.506 seconds, a 1.71x speedup (41.4% less wall time). Final objectives
agreed to about 1e-12, both exact final gradient norms were below 1e-4, and both
processes peaked at 18.078 MiB RSS.

Three-process threshold sweeps refined the policy. Exact-only median wall time
was 31.096 seconds. Threshold 0.03 converged in all runs with a 13.444-second
median, giving a 2.31x speedup and 56.8% wall-time reduction. Thresholds 0.1
and 0.3 both switched at iteration 52 and had roughly 18.45-second medians.
Threshold 0.03 switched at iteration 64, used 67 approximate and 8 exact
evaluations, and is the current catch-at-age recommendation.
