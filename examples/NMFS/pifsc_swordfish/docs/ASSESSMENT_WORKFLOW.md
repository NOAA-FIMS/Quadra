# Assessment workflows and comparison history

Run commands from the `examples/NMFS/pifsc_swordfish` case directory unless
another location is specified. This guide records the original comparison and
audit workflows; use [native input](NATIVE_INPUT.md) for the production C++ CLI.

## Run the model foundation and checks

From this directory:

```bash
make test
```

The audit/preparation commands below require Python 3 and a C++17 compiler.
The production native workflow requires only the built executable.
Preparation imports original observations and
controls, compiles the seasonal kernel, checks its conditional replay against
the supplied report, and checks Quadra AD derivatives and parser safeguards.
`make prepare` generates inputs without compiling; `make audit` runs the
historical conditional comparison, `make audit-biology` checks biological
predictions, and `make audit-population` checks the independently generated
historical trajectory at SS3 fitted parameter values. Generated inputs,
comparison CSVs and fitted outputs are written under `build/`.

See [MODEL_AND_COMPARISON.md](MODEL_AND_COMPARISON.md) for the schema, equations,
validation scope, fixed male-growth discrepancy and explicit SS3 3.30.20
spawning-timing convention. The replay calculates F-at-age, weights and mature weights from parameters.
The older conditional replay still uses reported apical F, recruitment and the
initial population. The population audit calculates all three from original
catches and the explicitly supplied fitted parameter values, including 57
recruitment deviations. These audits check predictions at supplied values; the
separate fitting executable estimates parameters.

## Evaluate and fit

```bash
make audit-likelihood
make objective-gradient
make fit-swordfish
```

The fit starts from the explicitly mapped SS3 estimates and optimizes 84
independent bounded parameters, with main recruitment deviations constrained to
sum zero and survey catchability profiled analytically. Recruitment deviations
are optimized jointly; this fit does not use Laplace marginalization.

The first fit converged at objective **289.299639305825**, gradient norm
**8.94e-5** (requested tolerance 1e-4). The supplied SS3 report prints **289.300**.
All 84 automatic derivatives passed centered numerical checks. This establishes
a common-start reproduction fit; multiple starts, uncertainty, reference points
and projections remain to be completed.

Results are written to `build/quadra_fit_summary.json`,
`quadra_fitted_parameters.csv`, `quadra_fitted_annual.csv`,
`quadra_fitted_seasonal.csv`, `quadra_fitted_index_fits.csv`,
`quadra_fitted_length_fits.csv`, and `quadra_fitted_likelihood.csv`.
New fits also write `quadra_fit_trace.csv`. `make fitted-report` regenerates the
fitted outputs from the saved parameter CSV without re-optimizing.
`make compare-fit` also checks all 47 annual SSB, recruitment, biomass and F
values against SS3 and writes `quadra_fit_comparison.csv` and `.json`.

## Reproduce the intake

From the repository root, using Python 3 (standard library only):

```bash
python3 examples/NMFS/pifsc_swordfish/scripts/import_ss3.py --source /path/to/supplied/files
```

The source directory must contain the eight filenames listed in the script.
Original files are copied byte-for-byte to `data/raw/`, which is excluded from
git. `reference/source_manifest.json` records their sizes and SHA-256 hashes.
The importer refuses to replace a differing raw file. It supports this supplied
SS3 3.30.20.00 report layout, not arbitrary SS3 versions.

`data/normalized/` contains report-derived CSVs with SS3 labels, missing tokens
(`_`, `NA`), original precision, and no unit conversions. Blank trailing fields
in derived quantities and parameters were absent in the original report.
These tables are benchmarks, **not independent raw observations for model fitting**.
The separate model-input parser imports catch, CPUE, length frequencies and
likelihood controls directly from `.dat` and `.ctl`.

| File | Contents |
|---|---|
| `fleets.csv` | 27 fleet definitions with original catch-unit codes |
| `parameters.csv` | 132 fixed and estimated parameters, phases, uncertainty, gradients |
| `derived_quantities.csv` | 261 quantities including SSB, recruitment, SPR, F ratios and reference points |
| `seasonal_time_series.csv` | 200 rows including virgin, initialization, historical and forecast eras |
| `annual_stock_status.csv` | 48 annual B/Bmsy and F/Fmsy rows; filter 1975–2021 for historical comparisons |
| `catch_fits.csv` | 3,648 fleet/season catch predictions and observations, including initialization |
| `index_fits.csv` | 164 index fit rows with observed/predicted indices, SE, likelihood and use flags |
| `length_fit_summary.csv` | 351 composition fit summaries; not bin-level frequencies |
| `likelihood_components.csv` | Total and component objective values, retaining SS3 conventions |

The script checks table widths, historical year/season coverage and unique
keys, and agreement of terminal SSB/SSB_MSY across report tables. It does not
validate biological equivalence, uncertainty coverage, or convergence beyond
recording the reference output. `reference/intake.json` records dimensions and
row counts.

## Confirmed input structure

- 1975–2021, four three-month seasons, two subseasons.
- Two sexes, ages 0 through 15+ (16 bins), one model area.
- 19 catch fleets and eight survey/index fleets; fleet names containing areas
  do not imply multiple population areas.
- Spawning month 5; recruitment settlement month 7, age 0, as specified by
  the control file's settlement assignment. Preserve SS3 age/timing semantics.
- Sex-specific age-dependent M and von Bertalanffy growth; length-based maturity
  and selectivity; fleet/time-block-specific selectivity and mirroring.
- Beverton–Holt recruitment, fixed h = 0.9 and sigmaR = 0.42, with early/main
  deviations and a fitted bias-adjustment ramp.
- Hybrid F method 3; catch observations mix biomass and numbers. Preserve fleet
  units and establish absolute scaling before combining removals.

## Reference files and interpretation

The supplied `Report.sso` identifies a November 12, 2026 rerun of the named
2023 assessment inputs using SS3 3.30.20.00 (September 30, 2022 build, ADMB 13.0).
It reports 86 active parameters, objective 289.3, and final gradient 5.81133e-05.
This is the supplied numerical benchmark; equivalence to the original 2023 run
has not been independently established.

| Quantity | Supplied Report.sso value |
|---|---:|
| SSB_2021 | 35,778.2 |
| SSB_unfished | 99,807.8 |
| SSB_MSY | 16,387.8 |
| annF_MSY | 0.178078 |
| Dead_Catch_MSY | 14,924.5 |

The official assessment PDF, Table S2 (printed page 5; PDF page 7), labels
95,732 as the average five-year dynamic unfished SSB baseline. This is distinct
from `SSB_unfished`; do not use one as the denominator for the other. Narrative
stock-status wording has internal inconsistencies, so compare labeled tables
and definitions rather than copying rounded prose ratios.

The summary document describes ten-year, five-scenario projections. The supplied
`forecast.ss` has forecast mode 0 and zero forecast years. Published projection
scenarios therefore require additional configuration. Keep any such setup
separate from the original input files. Also audit raw SPR versus 1-SPR:
starter and report descriptive labels are not entirely consistent.

`starter.ss` selects initialization from `ss.par`; that file was not supplied.
An exact reproduction of the original initialization requires it and the SS3
executable/build. Benchmark extraction and Quadra development can proceed.

## Implementation sequence

1. Implemented: import observations, sample sizes, weights, control parameters,
   blocks, mirrors, timing and unit definitions into a swordfish data schema.
2. Implemented: seasonal sex/age dynamics, length probabilities, plus-group
   growth, mature weights and selectivity, checked against SS3.
3. Implemented: initialization, recruitment bias adjustment, hybrid catches,
   CPUE and composition likelihoods, recruitment penalties, soft bounds, and
   objective matching at the supplied SS3 parameters.
4. First Quadra fit converged from the SS3 parameter estimates. Local curvature,
   parameter and selected derived uncertainty, and perturbed-start checks are
   described in [FIT_DIAGNOSTICS.md](FIT_DIAGNOSTICS.md). Recruitment
   treated as Laplace-integrated random effects changes the SS3 estimation
   target and must be labeled as a separate comparison.
5. Equilibrium MSY, SPR and biomass-target reference points match SS3. Dynamic
   unfished history and historical SPR also match. Five deterministic projection
   scenarios, source conflicts, and three interpretation sensitivities are in
   [PROJECTIONS.md](PROJECTIONS.md). Published bootstrap projections require
   their original forecast settings and bootstrap model outputs for exact replay. Local
   parameter intervals for all five scenarios are in
   [PROJECTION_UNCERTAINTY.md](PROJECTION_UNCERTAINTY.md).

The existing `pifsc_spatial_tuna` implementation supplies useful Quadra workflow
patterns, but its current population schema has no explicit sex dimension,
starts displayed ages at 1, and fits age compositions rather than this model's
length compositions. Reusing it unchanged would not reproduce the SS3 model.

## Operational release and code guide

Use `python3 scripts/run_release.py --output runs/your-new-run-name` for a clean,
isolated acceptance run and non-overwriting result publication. Low-level make
targets write scratch outputs. See [RELEASE_READINESS.md](RELEASE_READINESS.md)
for failure handling, covariance provenance, CI scope and remaining scientific
review requirements.

The [PDF operations and complete API guide](../output/pdf/swordfish_operations_and_api_guide.pdf)
and [searchable API reference](API_REFERENCE.md) document the swordfish
functions/classes and the shared Quadra changes. Regeneration instructions are
in [docs/README.md](README.md).
