# Reproducible native swordfish release

From this case directory, with the existing `.par`-start native JSON:

```bash
python3 scripts/run_native_release.py \
  --input build/native/swordfish-ss-par-start.json \
  --pdf-python /private/tmp/quadra_swordfish_docs_env/bin/python
```

That PDF environment is available in the current workspace. For a portable
setup, create a report environment under ignored `build/` and install ReportLab:

```bash
python3 -m venv build/report-env
build/report-env/bin/python -m pip install reportlab
make native-release PDF_PYTHON=build/report-env/bin/python
```

`NATIVE_INPUT` can select a different file path containing the same preserved
`.par` starting vector. `--output` selects a new bundle directory; otherwise the
runner creates a timestamped directory under ignored `runs/`. Existing destinations
are rejected. Progress appears in the terminal and each stage has its own log.

## What acceptance covers

1. Verify that all 74 native initial values and 57 deviations match the preserved
   `.par` import, before Quadra optimization.
2. Freeze framework/model source and a small set of comparison targets into a
   new workspace. Copy no raw SS3 files, generated observation factories,
   prebuilt executables, old fitted outputs or old covariance.
3. Compile the native executable and portable C++ tests from scratch. Run native
   schema/release safeguards, the optimizer-scaling chain-rule checks and five
   framework regressions.
4. Fit the baseline native model, independently check its original-coordinate
   gradient and calculate fresh two-step curvature/reference-point uncertainty.
5. Generate three deterministic perturbations of the `.par` starting vector,
   using seeds 101, 202 and 303 at 0.10, 0.25 and 0.50 local-SE scales along
   full baseline covariance directions. Preserve fixed settings and reconstruct
   the dependent final recruitment deviation.
6. Fit each perturbed input using that fresh covariance solely as an optimizer
   metric. Check convergence in original coordinates, objective recovery within
   `1e-5` and every independent coordinate within `0.01` local SE of the baseline.
7. Compare the baseline objective with the recorded SS3 objective and canonical
   Quadra fit, all ten reference points with printed SS3 targets, and all ten
   corresponding reference-point SEs with printed SS3 SEs.
8. Generate a two-page PDF summary, verify source/target identity, hash all
   artifacts and atomically rename the staging directory into the new destination.

This is native-path reproducibility and **local** start recovery. Jon's `.par`
file already contains a saved SS3 fitted vector. SS3 is not rerun, and these
perturbations do not prove global optimality. Missing SS3 bootstrap/projection
inputs and separate striped marlin work remain outside this acceptance scope.

## Bundle layout and provenance

| Path | Contents |
| --- | --- |
| `release_manifest.json` | Passed status, environment, stage timings, exact source and artifact hashes |
| `source_snapshot.json` | Frozen source and comparison-target hashes |
| `acceptance.json` | Baseline, numerical checks and local-start recovery results |
| `reference_point_comparison.json` | Ten SS3 point-estimate comparisons |
| `reference_SE_comparison.json` | Ten SS3 SE comparisons |
| `inputs/` | Exact baseline JSON and all seeded perturbation documents |
| `results/baseline/` | Fit, gradients, fresh covariance, reference points and uncertainty |
| `results/start-01/` through `start-03/` | Independently checked local fits and optimizer traces |
| `logs/` | Full build, test and model-stage output |
| `workspace/` | Model and bundled framework source, plus clean-build executables |
| `output/pdf/native_release_summary.pdf` | Actual accepted results and intervals |

Every fit records its exact input and executable hashes. Perturbed fits also
record the fresh optimizer-metric hash. The scaled optimizer trace uses scaled
gradient norms; `objective_gradient.csv` and final summaries always report the
original bounded-logit coordinates and independently checked gradients.
Uncertainty is calculated only for the baseline; local starts reuse its matrix
for optimization, not as newly estimated uncertainty.

Failed runs retain logs and partial artifacts in a clearly named private
`.failed-...` directory. Stage deadlines terminate the process group and no
failed bundle is published as passed. The final local bundle contains original
observations in normalized JSON and should be shared deliberately; it remains
under ignored `runs/` by default.

To reproduce from a bundle later, run its frozen `scripts/run_native_release.py`
with `--input` pointing to that bundle's `inputs/model.json`, a new `--output`
and a Python environment containing ReportLab. The frozen source snapshot contains
the framework and selected targets needed for another clean native run.

## Optional optimizer scaling

The native frontend additionally accepts:

```bash
./swordfish --input MY_PERTURBED_INPUT.json --output NEW_RESULTS \
  --optimizer-metric PRIOR_BASELINE_RESULTS/reference_coordinate_covariance.csv
```

The metric must be positive definite, labeled in exactly the model's coordinate
order and unchanged during the run. It changes optimizer steps through a constant
linear transformation, adding no statistical Jacobian or likelihood term.
Scaled convergence is followed by an independent original-gradient check.
The release runner uses only the covariance produced inside its own baseline run.

## Accepted local run

`runs/native-release-20261115-001` completed all ten stages from a clean
snapshot. The baseline objective is `289.299639305649`. The three perturbed
fits recovered `289.299639305628`, with final original-coordinate gradient
norms between `6.7e-8` and `9.2e-7`. Every coordinate returned within
`3.9e-6` local SE of the baseline. All ten SS3 point-estimate comparisons,
all ten SS3 SE comparisons and both numerical uncertainty steps pass.

All 2,275 artifact hashes and 2,183 frozen source/target hashes were verified
after publication. The two-page PDF was rendered and visually checked.
The [compact acceptance record](reference/native_release/acceptance.json),
[bundle summary](reference/native_release/summary.json) and
[PDF summary](output/pdf/native_release_summary.pdf) are available separately;
the complete observation-bearing bundle remains local under ignored `runs/`.
