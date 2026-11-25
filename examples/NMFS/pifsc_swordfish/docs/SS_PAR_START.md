# Initialize Quadra from the supplied SS3 parameter file

Jon's `ss.par` is an ADMB/SS3 parameter snapshot that can be read as a starting
vector. Its header also records an already fitted state:

```text
Number of parameters = 86
Objective function value = 289.299635060701
Maximum gradient component = 5.82152741601093e-05
```

Using this snapshot checks agreement when Quadra starts from that saved SS3
vector. The file does not establish the original pre-optimization assumptions,
and this experiment does not rerun SS3. Confirming recovery from a common
pre-fit vector would require that earlier vector and matching SS3 run evidence.
The separate [control-start experiment](CONTROL_START.md) starts from the
original `.ctl` values and zero recruitment deviations.

## Results from this file

| Quantity | Value |
| --- | ---: |
| SS3 objective recorded in `ss.par` | 289.299635060701 |
| Quadra final objective | 289.299639305649 |
| Quadra minus SS3 objective | 0.00000424495 |
| Quadra final gradient norm | 0.0000938583 |
| Difference from existing Quadra optimum | -0.000000000176 |

Both objective comparisons pass a declared absolute tolerance of `1e-4`.
The objectives agree to that tolerance; they are not bit-for-bit identical.
The mapped initial vector is already at essentially the same objective, but
Quadra still takes optimization steps to satisfy its own gradient criterion.
SS3's reported maximum gradient and Quadra's gradient use different coordinate
representations and should not be compared directly.

The largest final parameter difference from the saved Quadra reference is
`1.38e-5` reference standard errors. The largest annual spawning-biomass
difference is `0.1084` tonnes. These reference outputs and standard errors are
read only after fitting, for comparison.

## Mapping and supported scope

The strict parser maps `MGparm[1..20]`, `SR_parm[1..5]`, `init_F[1]`,
`Q_parm[1..8]` and `selparm[1..40]` to control-file structural order, including
expanded selectivity time blocks. It maps 20 early deviations to 1965-1984 and
37 main deviations to 1985-2021. This yields 74 native values and 57 deviations.
Unknown, duplicate, missing, incorrectly sized and nonfinite groups fail.
This importer supports this assessment layout, rather than arbitrary SS3 files.

The file supplies the fixed male growth value `83.4` directly; the control-file
value is `83.2`. No manual growth override is needed for the `.par` experiment.
All fixed values are taken from the supplied parameter snapshot.

The SS3 header's count of 86 is preserved as metadata. Quadra optimizes 84
independent coordinates, profiles Q analytically, and reconstructs the last
main recruitment deviation as minus the sum of the others. The input main
deviations sum to approximately `4.9e-13`; any reconstruction adjustment is
at that rounding scale. The Q values are imported, but Quadra's objective
profiles them at every evaluation. Dummy, zero forecast recruitment and
checksum entries are recorded separately; the checksum is not independently
recomputed. Forecast entries do not enter this historical fit.

## Run locally

From this case directory:

```bash
python3 scripts/run_ss_par.py /Users/matthew.supernaw/Downloads/ss.par
```

The original file has also been preserved locally at ignored
`data/raw/ss.par`. With that file present:

```bash
make fit-ss-par
```

Requires Python 3, a C++17 compiler (`CXX` selects it), repository headers,
original `.dat`/`.ctl` files under `data/raw/`, and saved reference summaries for
the post-fit comparisons. Neither preparation nor optimization reads
`Report.sso` or a saved covariance. The runner compiles the existing assessment
entry point against a factory generated only from `.dat`/`.ctl`, then supplies
native/deviation fixtures derived solely from `.par`.

Every run creates a new timestamped directory under `build/ss_par/`. It shows
optimizer progress in the terminal and saves `fit.log`, `initial_gradient.log`,
`initial_parameters.csv`, `input_manifest.json`, `comparison_summary.json`,
comparison CSVs, and fitted outputs under its own `build/` subdirectory. Existing
output directories are rejected. The canonical fit is not overwritten.

Reference points are now automatic fitted-model outputs in
`build/management_reference_points.csv`. Each `.par` run also writes
`benchmark_comparison.json` and requires all ten reference-point checks against
the saved printed SS3 targets to pass, in addition to objective convergence.
These comparisons happen after the model independently calculates its points.
The model API exposes the same calculation as `model.reference_points(x)`.
The standard historical fit, control-start fits and report regeneration also
write the table; no separate benchmark command is required.

The [saved experiment](../reference/ss_par/comparison_summary.json) and its
[hash manifest](../reference/ss_par/manifest.json) preserve the derived results;
the original `.par` file remains excluded from version control. No new
uncertainty or projection results are certified by this test. The subsequent
[integrated reference-point run](../reference/model_reference_points/benchmark_comparison.json)
passes all ten SS3 checks, with maximum relative error about `5.69e-6`.

The [one-page PDF report](../output/pdf/swordfish_ss_par_start.pdf) presents the
comparison. Regenerate it with `python scripts/build_ss_par_report.py` in the
documentation environment containing ReportLab.
