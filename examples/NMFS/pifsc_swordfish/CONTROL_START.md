# Starting from original SS3 inputs

Quadra can fit the swordfish assessment from the original `.dat` observations
and `.ctl` initial values. The control-only executable does not load
`Report.sso`, fitted parameter CSVs or a saved covariance. The existing canonical
fit still starts from SS3 estimates and remains unchanged.

## Results

| Configuration | Fixed male L1 | Joint objective | Full logit gradient norm |
| --- | ---: | ---: | ---: |
| Strict control inputs | 83.2 | 289.236730847211 | 7.23e-7 |
| Control starts, matched fixed growth | 83.4 | 289.299639305628 | 1.15e-8 |
| Existing SS3-initialized Quadra fit | 83.4 | 289.299639305825 | 8.94e-5 |

The control file fixes male `L_at_Amin_Mal_GP_1` at 83.2, while the supplied
fitted SS3 report uses 83.4. Strict and matched configurations therefore define
different objectives. Matching this one fixed setting permits a direct check
of initialization without borrowing any estimated free parameters.

The matched configuration recovers the existing optimum: objective difference
-1.97e-10, maximum parameter difference 1.39e-5 reference standard errors,
and maximum annual spawning-biomass difference 0.05893 tonnes. The saved
[comparison JSON](reference/control_start/comparison/comparison_summary.json)
contains the acceptance thresholds and input hashes. Reference standard errors
are read only after fitting to express parameter differences; they do not
provide optimizer scaling.

## Initialization and numerical method

1. Read all 74 native values from the original control initial-value column.
   Start all 57 historical recruitment deviations at zero. There are 84
   independent estimated coordinates after the main-deviation sum constraint.
2. Evaluate the literal control point. Its full likelihood fails with
   `Empty composition`: the initial recruitment scale cannot sustain the observed
   catches and leaves an empty predicted composition. This failure is recorded.
3. Adjust recruitment scale with an explicit catch-only feasibility objective.
   This starts from the literal controls; it does not alter the final likelihood
   or silently insert fitted values.
4. Activate coordinates using the control phases: 1, 8, 64 and 84 coordinates.
   These stages use the final weighted joint objective, rather than replaying
   SS3's phase-dependent likelihood weights.
5. Refine using linear scaling learned from centered differences of exact AD
   gradients at the cold checkpoint. Absolute eigenvalues and an eigenvalue
   floor make an optimization metric; this is not an inferential covariance.
6. Restore numerical sensitivity if a bounded logit has magnitude above 16,
   moving its starting probability to 0.0001 or 0.9999. The likelihood and
   bounds are unchanged. A final inward probe rejects apparent convergence if
   a saturated coordinate still has a decreasing inward direction.

The first small-gradient results pinned `SizeSel_fleet_8_p2` near its lower
bound of -10 and missed the optimum. The inward restart recovered the interior
solution (the matched reference is about -7.69). This is why gradient size alone
was insufficient for accepting those earlier results.

## Run and inspect

```bash
make test-control-start
make fit-control-start
```

The dedicated factory generator and executable require only the two original
`.dat`/`.ctl` inputs. An isolated build/test with those two files passed; the
saved [test log](reference/control_start/input_only_build_test.log) records it.
The comparison stage additionally reads saved reference summaries and outputs.
Python 3, a C++17 compiler and the repository's Eigen headers are required.

The runner creates a new directory under `build/control_start/`, records input
and binary hashes, runs both configurations, then compares their results.
Each fit automatically writes `management_reference_points.csv` using its own
fitted parameters; strict controls and matched fixed growth retain their distinct
configurations. This new output does not retroactively add reference-point
validation to the historical checkpoint snapshots below.
Optimizer progress appears in the terminal and is also saved in `strict.log`
and `matched_fixed.log` in that directory.
Existing output directories are rejected. `optimizer_iterations` records the
optimizer's iteration counter; line-search evaluation counts can be larger.
Historical snapshots retain the earlier counter label `evaluations`.

The executable also accepts its own saved cold checkpoint directory:

```bash
./build/swordfish_control_start matched-fixed NEW_OUTPUT 200 COLD_CHECKPOINT_DIRECTORY
```

It checks the original seed and fixed settings, then uses the latest recorded
refinement checkpoint, or the phase-4 checkpoint. This is a local experiment
resume interface, not a general SS3-output importer. Keep all ancestor outputs
when using it.

## Evidence and limits

The accepted experiment resumed its own phase-4 checkpoints (`003`), refined
them (`004`), then recovered the saturated coordinate (`005`). All three stages
are preserved under [reference/control_start](reference/control_start/manifest.json),
including initial values, stage traces, metric diagnostics and final predictions.
The manifest hashes saved evidence. Earlier stages ran earlier revisions of the
new executable, so this is not an immutable release snapshot or a fresh
end-to-end execution of the final source revision. `make fit-control-start`
provides that fresh-run path; the resumed results do not certify its runtime or
identical trajectory.

The existing warm-start covariance and release artifacts have not been
regenerated for these source changes. No cold-fit uncertainty or projections
are certified here. Missing SS3 bootstrap/projection inputs and the separate
striped marlin implementation remain outside this initialization experiment.

The [PDF report](output/pdf/swordfish_control_start.pdf) presents the results and
annual biomass comparison. Regenerate it with the documentation Python environment:

```bash
python scripts/build_control_start_report.py --evidence reference/control_start/comparison
```
