# Native input preparation

## North Pacific swordfish

The production converter is part of the C++ `pifsc_swordfish` executable.
From the repository root, build once and convert:

```bash
make -C examples/NMFS/pifsc_swordfish pifsc_swordfish
examples/NMFS/pifsc_swordfish/build/pifsc_swordfish --convert-ss3 \
  --data examples/NMFS/pifsc_swordfish/data/raw/swo2023_v004.dat \
  --control examples/NMFS/pifsc_swordfish/data/raw/swo2023_v007.ctl \
  --parameters examples/NMFS/pifsc_swordfish/data/raw/ss.par \
  --output examples/NMFS/pifsc_swordfish/build/native/imported-swordfish.json
```

Only the supported North Pacific swordfish SS3 3.30.20 layout is accepted.
Source filenames and directories may vary. Runtime needs only the built binary;
no Python, `Report.sso`, forecast file or covariance is required.

The optional `--parameters` supplies initial values and historical deviations
from `.par`, including fixed biology, and records saved fit metadata in provenance.
Omitting it keeps literal control starts and zero deviations; those starts require
the separate catch-feasibility workflow before fitting.

Conversion writes self-contained JSON and a `.provenance.json` sidecar with
source hashes. Existing outputs are protected. Then evaluate or fit directly:

```bash
examples/NMFS/pifsc_swordfish/build/pifsc_swordfish \
  --input examples/NMFS/pifsc_swordfish/build/native/imported-swordfish.json \
  --output examples/NMFS/pifsc_swordfish/build/native/imported-evaluation \
  --mode evaluate
```

Omit `--mode evaluate` to fit; add `--uncertainty` for reference-point intervals.
Use new output paths for each run. See the case's
[native input guide](../examples/NMFS/pifsc_swordfish/docs/NATIVE_INPUT.md).

The former Python launcher in this directory has been removed. Independent
Python conversion logic lives in the case's `tests/support/` for parity testing;
research, SS3 comparison and release/reporting tools remain under `scripts/`.
