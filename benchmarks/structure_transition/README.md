# Quadra structure-transition benchmark

This Quadra-only diagnostic holds the latent Gaussian workload constant while
moving from tridiagonal through wider bands to an O(n)-edge irregular sparse
precision matrix. It records automatic structure/backend selection, warm
persistent numeric-update runtime, log-determinant parity, and isolated process peak
RSS. Timing uses the requested repetitions; RSS is measured separately with
one warm persistent update so allocator high-water history does not distort the
structural memory comparison.
The CSV also reports `repeated_peak_rss_mib` for the complete timed loop, which
guards against allocator growth across persistent updates.

The compressed Hessian pattern is retained between timed evaluations; only its
numeric values are changed. This prevents input-matrix allocation from being
mistaken for backend factorization cost.

Run from the Quadra repository root:

```bash
./benchmarks/structure_transition/run_benchmark.sh 50 250,1000,5000
```

The irregular case retains a first-order Markov backbone and adds deterministic
long-range chords. It is sparse but has large bandwidth, so it exercises the
generic sparse fallback instead of a specialized banded factorization.

## Current result

At 5,000 latent effects:

| structure | backend | nnz | ms/update | peak RSS MiB |
|---|---|---:|---:|---:|
| tridiagonal | tridiagonal | 14,998 | 0.031 | 2.83 |
| banded-5 | banded | 54,970 | 0.204 | 6.20 |
| banded-10 | banded | 104,890 | 0.522 | 10.72 |
| banded-32 | banded | 323,944 | 3.306 | 31.11 |
| irregular sparse | sparse LDLT | 16,428 | 2.011 | 4.02 |

All automatic backend choices were correct and all log determinants matched an
independent sparse LDLT reference within `3e-12`.

The persistent sparse backend now uses the shared LDLT factorization cache. In
the 50-update run, every irregular case performs one symbolic analysis and 51
numeric factorizations (the initial evaluation plus 50 updates). At 1,000
effects this reduced the irregular update from about 0.136 ms to 0.050 ms; at
5,000 effects it reduced the update from about 2.57 ms to 2.01 ms.
