# Run the model directly from normalized input

From this case directory:

```bash
make native-swordfish
./swordfish --input build/native/swordfish.json --output build/native/my-fit
```

The local example `build/native/swordfish.json` is ready to use. It contains
original observations/settings and explicit initial values from the existing
Quadra fit. It is a self-contained input document, rather than a reference to
SS3 files. The output directory must be new.

The compiled backend is reusable. Changing the JSON input does not require
recompilation. Runtime needs Python 3 and the built executable; building needs
a C++17 compiler and this repository's headers. The native runner does not
read `.dat`, `.ctl`, `.par`, `Report.sso`, saved reference outputs or covariance.
SS3 comparison tools remain optional validation utilities.

To evaluate an initial model state without fitting:

```bash
./swordfish --input build/native/swordfish.json \
  --output build/native/my-evaluation --mode evaluate
```

Both modes write annual/seasonal predictions, observation fits, likelihood
components, native parameter values, exact AD gradients and management reference
points. Fit mode also saves the optimizer trace. Progress appears in the terminal.
Each successful run preserves the exact JSON as `model_input.json` and records
input/binary SHA-256 hashes in `run_status.json`. Evaluation mode reports
`converged: false` because no optimization was requested; its run status is
`completed` when evaluation succeeds. Failed fits return a nonzero exit status
and do not produce final fitted reports.

Options include `--max-iterations 200` and `--gradient-tolerance 0.0001`.
Add `--uncertainty` in fit mode for freshly calculated reference-point SEs,
CVs and pointwise 95% intervals. See
[REFERENCE_POINT_UNCERTAINTY.md](REFERENCE_POINT_UNCERTAINTY.md) for the full
correlated delta method, numerical acceptance checks and output files.
The default is the existing joint likelihood, including recruitment penalties;
no Laplace marginalization or new reference-point uncertainty is added.

## Native schema

The root object contains exactly:

| Field | Meaning |
| --- | --- |
| `schema` | `quadra.swordfish.v1` |
| `units` | Explicit numbers, weight, biomass, length and mortality units |
| `model` | Calendar, biology, fleets, observations, selectivity and likelihood settings |
| `parameters` | Named native parameter values, bounds, estimation phases and block metadata |
| `recruitment_deviations` | Explicit year/value records covering all historical deviations |

Parameters use named fields: `name`, `lower`, `upper`, `value`, `phase`,
`prior_type`, `block`, `block_function`. Positive parameter phases mean estimated;
nonpositive phases mean fixed. Only the currently supported no-prior setting
(`prior_type: 0`) is accepted. Selectivity block functions retain their model
codes: 0 multiplicative on the exponential scale, 1 additive, 2 replacement.

Recruitment settings name the early/main periods, phase settings, deviation
bounds and bias-ramp years/fraction. The historical deviations are jointly
estimated by this model, with the final main deviation dependent on the sum
constraint. The input main deviations must already sum to zero within `1e-8`;
the final value is reconstructed at rounding precision during model decoding.

Biology and observations use named records. Fixed-width settings tables retain
the backend's existing code conventions; this version preserves the supported
swordfish model rather than providing an arbitrary-species schema. It supports
two sexes, ages 0-15+, four seasons, two subseasons, one area, 19 fishing fleets
and eight index fleets, and the supplied combined-sex composition layout.
`scripts/native_input.py` defines required fields, types and table widths.
Unknown/missing fields, duplicate keys, nonfinite numbers, malformed arrays,
invalid estimated starts, missing parameter names and inconsistent deviation
years/sums fail validation. C++ also validates the loaded model state.

The literal low-recruitment control starts can produce an empty predicted
composition. The generic native runner rejects an unevaluable initial state;
it does not silently substitute a fitted vector or run the specialized cold
feasibility procedure. Use the explicit [control-start workflow](CONTROL_START.md)
for that experiment. The ready-to-use native example deliberately uses the
existing fitted vector, so this test establishes input-path equivalence, not
new evidence of convergence from arbitrary starts.

## Optional one-time preparation

To recreate the example from the existing parsed-input JSON and a complete
native parameter CSV:

```bash
python3 scripts/export_native.py \
  --source-normalized build/swordfish_inputs.json \
  --initial-parameters reference/model_reference_points/quadra_fitted_parameters.csv \
  --output build/native/new-model.json
```

Alternatively, use `--ss3-directory data/raw` instead of `--source-normalized`.
The exporter is separate from the runner and is never imported at runtime.
Omitting `--initial-parameters` explicitly selects control initial values and
zero deviations, with the feasibility limitation above. Existing files are not
overwritten. Once exported, the native document is all the model input needed.

Keep this observation-bearing JSON local unless its publication is intended.
The generated example remains under ignored `build/`; it is not added to Git.
The older CSVs under `data/normalized/` are assessment output summaries, not
a complete replacement for observations and model configuration.

## Validation and evening stopping point

The native evaluation/fit recovered objective **289.299639305649**. Annual,
seasonal, index, composition, likelihood and reference-point output tables
agree with the previous fit to better than `6e-14` relative error (using a
denominator of at least one). An isolated copy ran without SS3 sources,
reference fixtures or generated input factories, without invoking a compiler. Changing R0 by
1% in that isolated JSON changed unfished SSB by 1% using the same binary.

The [saved validation](reference/native_input/validation.json) records these
checks. Native-format rejection tests are included in the Python suite. This
is a native entry-point validation, not a refreshed full release/covariance run.
The existing canonical fit and historical validation snapshots are preserved.
For the next session, start with the two commands at the top of this document.
