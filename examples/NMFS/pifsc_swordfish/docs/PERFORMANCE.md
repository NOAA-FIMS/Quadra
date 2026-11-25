# Fixed-biology evaluation workspace

The first performance change is limited to biological quantities whose
parameters are fixed. Scientific equations remain in the existing biology,
population and likelihood functions.

`FixedBiology` in `include/swordfish/fixed_biology.hpp` owns a small native
parameter snapshot and stores ordinary numeric length weights, mature weights,
time-specific age-length probabilities and beginning-season individual weights.
Population and observation calculations share this workspace during one
likelihood evaluation. It is recreated on every evaluation, so changed inputs
are reflected immediately and no AD variables survive across tape boundaries.

The dependency guard is conservative: estimating any biological parameter
selects the original differentiable path for all these quantities. Estimated
selectivity curves and population dynamics retain their derivatives. Exactly
zero numeric probabilities can skip zero contributions; AD-valued probabilities
retain the full calculation, including when their current value is zero.

## Readability and review

The workspace has four named operations: construction, `enabled()`,
`length_keys(time)` and `weights_at_age(time)`. The existing integration equations
remain visible in `seasonal_biology` and `sampled_lengths`, with one branch to
choose numeric versus differentiable inputs. No new command-line setting is
required. The model's `use_fixed_biology_cache` switch is retained for validation
and benchmarking; biological estimation flags control the fallback automatically.

## Measurement

A small benchmark on this machine alternates original/cache order in the same
`-O2` binary. Each sample averages ten objective evaluations, six AD-gradient
evaluations and ten reference-point calculations; the reported comparison uses
the median of two samples per path. It is an evaluation benchmark rather than
a guaranteed end-to-end fit speedup or a cross-machine result. Repeated timings
and measurement details are saved with the validation evidence.

The final local measurements were:

| Operation | Original path | Fixed-biology workspace | Time reduction |
| --- | ---: | ---: | ---: |
| Objective | 6.74 ms | 3.72 ms | 44.8% |
| AD gradient | 420 ms | 349 ms | 17.0% |

[Benchmark samples and method](../reference/performance/benchmark_comparison.json)
are preserved with the measurement harness in `benchmark_source.txt`.

Reference-point calculations are already much cheaper than AD gradients. This
change targets repeated objective/gradient work, including the 336 gradient
evaluations needed by the two-step uncertainty calculation. It preserves both
curvature/Jacobian step checks.

## Validation

`make test-fixed-biology` compares cached versus original objective values and
every gradient component at the initial vector and a perturbed vector. It also
changes a fixed biological value between calls to check input refresh, then
estimates female growth K and checks the automatic fallback and its numerical
derivative. The test is included in the main `make test` target.

A fresh native fit with reference-point uncertainty checks fitted outputs,
reference points and SEs against the earlier validated run. Curvature remains
unregularized and must pass the same two-step positive-definiteness/stability
checks. Earlier release/covariance snapshots remain historical evidence; this
change does not refresh a full immutable release bundle.

The [fresh uncertainty comparison](../reference/performance/cache_validation.json)
recovers the previous objective to about `1e-13`, reference estimates within
`2.5e-12` relative difference and reference SEs within `1.2e-5` relative
(about 0.0012%). All ten saved SS3 reference-point SE checks still pass.
The full model test suite, including all 84 AD gradient checks and all 17
Python tests, passes. Run commands and scientific tolerances are unchanged.
