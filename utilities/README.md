# Input preparation utilities

These commands prepare native Quadra inputs from external formats. Model
execution remains separate and needs only the resulting native input document.
The CLI launchers live here; reusable, assessment-specific parsers remain with
their case implementations.

## North Pacific swordfish

From the repository root:

```bash
python3 utilities/ss3_to_pifsc_swordfish.py \
  --data examples/NMFS/pifsc_swordfish/data/raw/swo2023_v004.dat \
  --control examples/NMFS/pifsc_swordfish/data/raw/swo2023_v007.ctl \
  --parameters examples/NMFS/pifsc_swordfish/data/raw/ss.par \
  --output examples/NMFS/pifsc_swordfish/build/native/imported-swordfish.json
```

The source files can have arbitrary names and live in different directories.
Only the supported NP swordfish SS3 3.30.20 assessment layout is accepted;
this command is not a general SS3 converter. Python 3 is sufficient: no compiler,
`Report.sso`, forecast file or covariance is required for conversion.

The optional `--parameters` imports all values and historical recruitment
deviations from `.par`, including fixed biology. Its saved fit metadata is
preserved as provenance; it does not establish an original pre-fit vector.
Omitting this option uses control initial values and zero deviations. The literal
control point for this assessment cannot evaluate the full likelihood without
the separate, documented catch-feasibility stage.

The output is self-contained native JSON plus an adjacent `.provenance.json`
sidecar with source hashes, initialization details and counts. Existing JSON or
sidecars are protected. No source input is changed. Observation-bearing outputs
should normally be written under ignored `build/` unless sharing is intended.

Then run from the swordfish case directory:

```bash
cd examples/NMFS/pifsc_swordfish
make native-swordfish
./swordfish --input build/native/imported-swordfish.json \
  --output build/native/imported-fit
```

Add `--uncertainty` for reference-point intervals. Choose a new output directory
for each run. Existing preparation scripts remain available, including
`scripts/export_native.py` and `scripts/run_ss_par.py`.
