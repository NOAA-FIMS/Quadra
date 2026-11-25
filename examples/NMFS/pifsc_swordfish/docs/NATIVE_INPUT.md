# Run the model directly from normalized input

The compiled executable shows its full help menu with `--help`, `-h`, or no
arguments. Help needs no input files:

```bash
./build/pifsc_swordfish --help
```

The menu covers SS3 conversion, direct JSON model runs, optimizer options,
examples, provenance and the legacy stdin interface. All runtime work runs in C++.

From this case directory:

```bash
make pifsc_swordfish
./build/pifsc_swordfish --input build/native/swordfish.json --output build/native/my-fit
```

The local example `build/native/swordfish.json` is ready to use. It contains
original observations/settings and explicit initial values from the existing
Quadra fit. It is a self-contained input document, rather than a reference to
SS3 files. The output directory must be new.

The compiled backend is reusable. Changing the JSON input does not require
recompilation. Runtime needs only the built executable; building needs a C++17
compiler and this repository's headers. JSON and SHA-256 dependencies are vendored
under `include/swordfish/third_party/`; no network access is needed to build.
`./swordfish` is a shell convenience launcher for the same executable. Model evaluation and fitting
from JSON do not read `.dat`, `.ctl`, `.par`, `Report.sso`, saved reference
outputs or covariance. The optional `--convert-ss3` mode performs separate
one-time intake. SS3 comparison tools remain optional validation utilities.

To evaluate an initial model state without fitting:

```bash
./build/pifsc_swordfish --input build/native/swordfish.json \
  --output build/native/my-evaluation --mode evaluate
```

Both modes write annual/seasonal predictions, observation fits, likelihood
components, native parameter values, exact AD gradients and management reference
points. Fit mode also saves the optimizer trace. Progress appears in the terminal.
Each run that creates an output directory preserves the exact JSON as
`model_input.json` and records
input/binary SHA-256 hashes in `run_status.json`. Evaluation mode reports
`converged: false` because no optimization was requested; its run status is
`completed` when evaluation succeeds. Failed fits return a nonzero exit status
and do not produce final fitted reports. Failed runs after output creation also
record `status: failed` and their error in `run_status.json`. Malformed inputs
fail before creating an output directory.

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
`include/swordfish/native_json.hpp` defines required fields, types and table
widths; `scripts/native_input.py` retains the equivalent Python validation API
for optional analysis utilities.
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

Convert the supported North Pacific swordfish SS3 3.30.20 inputs directly with
the swordfish command:

```bash
./build/pifsc_swordfish --convert-ss3 \
  --data data/raw/swo2023_v004.dat \
  --control data/raw/swo2023_v007.ctl \
  --parameters data/raw/ss.par \
  --output build/native/from-ss3.json
```

Run `make test-ss3-conversion` for the converter and command-line checks
(using the local observation/control fixtures).

Place `--convert-ss3` first. Use `--convert-ss3 --help` for the conversion flags.
Only the explicit `.dat`, `.ctl`, and optional `.par` are read; filenames may
be arbitrary. `Report.sso`, forecast files, reference tables, covariance, and
supporting PDF/Word documents are unnecessary for this conversion.

The converter is compiled into the executable. A relocated binary needs no
`scripts/` directory or Python installation. The former `SWORDFISH_PYTHON` and
`SWORDFISH_CONVERTER` environment variables are no longer used.

Conversion writes the native JSON and a `.provenance.json` sidecar with source
hashes and initialization metadata. Existing outputs or sidecars are protected.
Omit `--parameters` to preserve literal control starts and zero historical
deviations. Those starts retain the catch-feasibility limitation described
above; conversion does not adjust starts or run a fit. Once conversion finishes,
run directly from the new JSON:

```bash
./build/pifsc_swordfish --input build/native/from-ss3.json \
  --output build/native/from-ss3-fit
```

SS3 intake is invoked only by the explicit conversion option. Normal JSON
model runs do not load the SS3 converter. See the [utilities guide](../../../../utilities/README.md)
for commands from the repository root.

## Python tools retained for validation and research

`pifsc_swordfish` is the production entry point. The duplicate Python JSON runner,
exporter CLI and converter launchers have been removed.

An independent converter remains under `tests/support/ss3_oracle.py`, with its
schema mapping in `tests/support/native_export.py`. These modules have no CLI;
the tests use them to compare C++ conversion against a separate implementation.
`scripts/native_input.py` remains shared by tests and release tooling to validate
and perturb research inputs. Original SS3 parsers also support the comparison
and control-start experiments. Python reporting and release tools remain optional.

Keep this observation-bearing JSON local unless its publication is intended.
The generated example remains under ignored `build/`; it is not added to Git.
The older CSVs under `data/normalized/` are assessment output summaries, not
a complete replacement for observations and model configuration.

## Validation

The native evaluation/fit recovered objective **289.299639305649**. Annual,
seasonal, index, composition, likelihood and reference-point output tables
agree with the previous fit to better than `6e-14` relative error (using a
denominator of at least one). An isolated copy ran without SS3 sources,
reference fixtures or generated input factories, without invoking a compiler. Changing R0 by
1% in that isolated JSON changed unfished SSB by 1% using the same binary.

The [saved validation](../reference/native_input/validation.json) records these
checks. Native-format rejection tests are included in the Python suite. This
is a native entry-point validation, not a refreshed full release/covariance run.
The existing canonical fit and historical validation snapshots are preserved.
For the next session, start with the two commands at the top of this document.

## Standalone C++ validation

Run `make test-native-cli` for portable C++ schema/SHA-256 checks and Python-driven
CLI regression tests. With the local SS3 fixtures, the tests compare every native
model field against the independent Python converter, check source/output hashes,
exercise malformed inputs and failed-fit status, and evaluate a relocated binary
with an empty PATH and no Python/scripts. Model JSON is semantically equivalent;
JSON formatting and insignificant summation rounding in metadata may differ.

The compiled command also retains the private positional stdin interface for
existing analysis utilities. New workflows should use named `--input` and
`--output` arguments. Specialized control-start experiments, SS3 comparison
audits and PDF release tooling remain separate optional Python utilities.
