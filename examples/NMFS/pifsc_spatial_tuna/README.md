# Quadra tuna assessment

For the complete scientific and computational specification—including every
state equation, parameter transform, likelihood, prior, simulation step,
Laplace calculation, sampler, reference point, projection, output, and known
limitation—read [MODEL_AND_QUADRA.md](MODEL_AND_QUADRA.md). This README is the
operator's quick start; that document is the auditable definition of the
assessment.

The assessment is configured at runtime with
`config/tuna_assessment.conf`. The format is dependency-free `key = value`
text; blank lines, `#` comments, and section labels are allowed.

This directory is a self-contained Quadra fisheries example. Run the commands
below from `examples/NMFS/pifsc_spatial_tuna`; its Makefile resolves the
Quadra root three directories above it.

The synthetic assessment does not require the exploratory OPAL `.rda` files
from the development workspace. Those data and the papers used to audit them
are intentionally not vendored into Quadra; `scripts/audit_opal_wcpo_bet.R`
remains as an optional importer for users who supply `data/opal_raw` locally.

```bash
make fit-advanced-tuna
make sample-advanced-tuna
```

Use a different configuration without recompiling:

```bash
make fit-advanced-tuna CONFIG=config/my_assessment.conf
```

The executable also accepts the file directly:

```bash
build/advanced_tuna_spatial_assessment_example \
  --config config/tuna_assessment.conf
```

Precedence is:

```text
compiled defaults < configuration file < QUADRA_TUNA_* environment variables
```

The resolved settings used by each run are written to
`build/assessment_outputs/data/effective_configuration.csv`. Machine-readable
assessment products are kept in `build/assessment_outputs/data`; the Markdown
report and figures are kept in `build/assessment_outputs/report`.

`sampling.gradient_workers` controls exact marginal-gradient parallelism
inside each chain. The default is one because chains themselves run in
parallel; higher values can oversubscribe the CPU and should be benchmarked on
the target machine.

`sampling.tape_rebuild_interval` bounds reuse of Quadra's persistent
structure-aware random-Hessian tape. The tuna workflow defaults to rebuilding
for every marginal evaluation (`1`) because longer reuse triggered stale
frozen-topology reads in long NUTS runs. The tape is still reused across the
inner Newton iterations within an evaluation.

Sampling performance counters are written to
`marginal_sampler_profile.csv` beside the other machine-readable assessment
outputs. They separate inner random-effect mode-solving time from exact
Laplace-gradient time and report the mean Newton iterations per evaluation.

`sampling.method` selects `pcn` (the clean-checkout default) or `nuts`. The exact pCN
Metropolis sampler works in the marginal-Hessian-whitened coordinates and
uses one posterior evaluation per transition. Its proposal scale is frozen
before retained sampling, so draws target the same Laplace-marginal posterior
as AD-NUTS.
For reproducible cross-chain geometry, the tuna configuration defaults to a
fixed `sampling.pcn_beta = 0.5`; experimental warmup adaptation can be enabled
with `sampling.pcn_adapt_beta = true`.
For pCN, marginal evaluations omit the exact fixed-effect gradient and mixed
derivatives. The validated whitening matrix is cached in
`marginal_whitening_cache.csv`; its fingerprint includes the assessment
inputs, fitted parameters, active parameter set, and whitening constants.

`sampling.method = covariance_rw` uses a frozen full-covariance random-walk
proposal learned from the configured exact pilot draws. It is symmetric, so
the retained kernel needs only the exact target-density Metropolis correction;
the tuna calibration uses `sampling.proposal_rw_scale = 0.2`.

Experimental frozen nonlinear proposals are available as `transport_gmm` and
`transport_kde`. Both use exact independence-Metropolis correction, but the
current 1,000-draw, 23-dimensional pilot set is too sparse for these proposals
to pass sampler health checks. They are not the configured production default.
`transport_poly` provides a dependency-free quadratic masked-autoregressive
flow with held-out validation, but it is likewise experimental; the current
pilot fit does not achieve usable exact-correction acceptance.

After a pCN sampling run has created pilot posterior draws, the neural transport
training stage is run with `make train-transport-flow`.
It invokes a standalone Eigen/C++ trainer that fits a deterministic four-seed
ensemble of eight-layer RealNVP affine coupling flows using whole-chain holdout
validation, exact deduplication, early stopping, explicit reverse propagation,
gradient clipping, AdamW, and serialized forward/inverse round-trip checks.
Frozen weights and normalization are written to four versioned 272-KB QFLOW
archives. Training and assessment use only C++, Eigen, and the standard
library; Python, nflows, TorchScript, and LibTorch are not used.

The native archive is loaded by `sampling.method = transport_flow` or
`transport_isir`. Its `<artifact>.manifest` is verified before weights are
loaded, including the manifest version, artifact-content hash, assessment
data/control fingerprint, fitted-mode and whitening-geometry fingerprint,
active dimension, and exact parameter names/order. The manifest preserves the
flow architecture, training sources and their hashes, training seed, held-out
validation likelihood, inverse round-trip error, and portable/reference parity
errors. Missing, modified, stale, or non-native artifacts fail closed before
sampling.

`sampling.method = transport_isir` uses the configured comma-separated
`sampling.transport_models` as an equal-weight frozen proposal mixture inside
an exact iterated sampling-importance-resampling kernel. Each transition keeps
the current state as one candidate, draws `sampling.transport_candidates - 1`
new candidates, evaluates exact target-to-proposal importance weights, and
resamples one candidate. This preserves the same Laplace-marginal posterior
while reducing the sticky-chain behavior of a single independence-Metropolis
proposal. Its exact target-evaluation cost per transition is recorded in
`sampler_identity.csv`.
Sampler health is written to the sampler-neutral `sampler_summary.csv`.
`nuts_summary.csv` remains as a temporary compatibility copy for downstream
consumers that have not migrated.

`make check-transport-gradient` compares selected analytical RealNVP gradients
with centered finite differences and checks serialized inverse round trips.
`make check-transport-reproducibility` trains the same small flow twice and
requires byte-identical QFLOW archives.

After marginal sampling, every retained fixed-effect draw receives a
conditional Laplace reconstruction of the recruitment random effects. The
workflow first solves their conditional mode and then draws from the local
Gaussian approximation using the random-effect Hessian as the precision
matrix. `posterior_random_effect_draws.csv` records both values for audit, and
`posterior_reconstruction_summary.csv` records failures. A reproducibly
thinned set of these joint draws, controlled by `sampling.management_draws`,
feeds `posterior_reference_points.csv` and
`posterior_projection_draws.csv`; no management calculation silently uses the
fixed-effect draw with latent effects left at the fitted mode.
`make check-posterior-assessment` enforces sampler convergence, complete and
finite latent reconstruction, valid reference points, four complete
biologically admissible projection scenarios, and separation of CSV data from
report artifacts. `make sample-advanced-tuna` runs this gate before building
the report.
