# Independent conversion oracle

`ss3_oracle.py` and `native_export.py` preserve the original Python conversion
logic as a reference for the C++ converter. They have no command-line entry point
and are not imported by production code. Tests add `scripts/` to their import
path for the shared original SS3 parsers, schema validator and hashing helpers.

Run from the case directory:

```bash
make test-native-cli
```

Conversion parity tests need the local SS3 fixtures in `data/raw/`; portable
C++ schema/SHA-256 and Python schema/release checks do not need those sources.
The release snapshot copies `tests/` recursively, including this support package.
