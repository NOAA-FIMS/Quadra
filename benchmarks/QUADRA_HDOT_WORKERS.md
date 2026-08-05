# Quadra Hdot worker study

This is a Quadra-only optimization benchmark. Its model uses the conventional
TMB Gaussian AR(1) state-space parameterization:

- unconstrained AR correlation;
- log process standard deviation;
- log observation standard deviation;
- process mean;
- a random latent state for every observation.

The covariance parameters affect the random-effect Hessian and are active Hdot
directions. The mean is intentionally Hessian-inactive, allowing Quadra's
automatic direction diagnostics to demonstrate pruning.

Run:

```sh
make benchmark-hdot-workers
```

The CSV output compares one worker, two workers, and automatic worker selection
at 100, 300, and 1,000 latent states. It reports setup and warm evaluation
times, exact-gradient phase timings, speedups, Hessian structure, backend and
solver recommendations, expected complexity, symbolic reuse, tape rebuilds,
and serial-versus-worker numerical differences.

The benchmark is deliberately not an RTMB timing comparison. TMB-compatible
model forms keep the study relevant to later parity work, while Quadra's
diagnostics identify optimization opportunities independently.

## Current diagnostic result

The worker tapes are initialized from one canonical recorded tape. Additional
workers reuse its discovered graph instead of re-recording the model. The
immutable flat-edge registry and per-vertex propagation-slot plan have shared
ownership; worker primals, local derivative weights, adjoints, tangents, and
flat value arrays remain independent. The `topology_owners` diagnostic verifies
that every worker references the same immutable topology.

Replay operation metadata is also separated from mutable vertices and shared:
opcode, operand and reverse-edge destination vertex IDs, and constants live in
one immutable operation table. Propagation uses transient edge views that join
this shared topology to worker-local derivative weights. `operation_owners`,
`operation_count`, `avoided_operation_bytes`, and
`avoided_edge_destination_bytes` make this visible in the CSV. At 1,000 states
the graph contains 24,018 operations; three workers avoid about 1.15 MB of
duplicated operation metadata, including about 0.38 MB of logically duplicated
edge destination IDs.

Active-direction discovery also runs on the canonical persistent tape, so its
work is retained rather than discarded.

For the 1,000-state model this reduced serial setup from roughly 355 ms to
67 ms. Three-worker setup is about 69 ms, only about 3% above serial, while the
warm exact-gradient evaluation improves from about 4.0 ms to 2.7 ms. Hdot
itself improves from about 2.9 ms to 1.4 ms. Exact timings vary with thread
scheduling, so the CSV includes both phase and end-to-end speedups.
