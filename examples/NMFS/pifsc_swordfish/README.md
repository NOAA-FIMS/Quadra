# PIFSC swordfish

A Quadra C++ implementation of the 2023 PIFSC/ISC North Pacific swordfish
assessment, with seasonal two-sex population dynamics, SS3 comparison checks,
likelihood fitting, management reference points and local uncertainty.

## Build and run

Run from this case directory. Building requires a C++17 compiler and the Quadra
headers in this repository. Conversion, evaluation and fitting need only the
built executable; Python is not required at runtime.

```bash
make pifsc_swordfish
./build/pifsc_swordfish --help

./build/pifsc_swordfish --convert-ss3 \
  --data data/raw/swo2023_v004.dat \
  --control data/raw/swo2023_v007.ctl \
  --parameters data/raw/ss.par \
  --output build/native/from-ss3.json

./build/pifsc_swordfish --input build/native/from-ss3.json \
  --output build/native/evaluation --mode evaluate
```

Omit `--mode evaluate` to fit; add `--uncertainty` for reference-point intervals.
Use new output paths for each run. The `.par` is optional; literal control starts
require the separate catch-feasibility procedure before fitting. Raw SS3 files
are local, ignored inputs. The `./swordfish` shell launcher forwards to the same
compiled executable.

## Documentation

Start with the [documentation index](docs/README.md).

- [Native input, conversion and CLI](docs/NATIVE_INPUT.md)
- [Reference-point uncertainty](docs/REFERENCE_POINT_UNCERTAINTY.md)
- [Reproducible native releases](docs/NATIVE_RELEASE.md)
- [Model equations and SS3 comparisons](docs/MODEL_AND_COMPARISON.md)
- [Control-start experiments](docs/CONTROL_START.md)
- [Assessment workflows and comparison history](docs/ASSESSMENT_WORKFLOW.md)
- [Searchable API reference](docs/API_REFERENCE.md)

## Checks

```bash
make test-native-cli
```

This runs portable C++ schema/hash checks and Python-driven CLI/parity tests.
Conversion checks use the local observation/control fixtures; optional research,
SS3 comparison and reporting tools use Python. `make test` runs the broader
assessment audit suite with the original SS3 inputs.
