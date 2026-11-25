# PIFSC swordfish documentation

Run commands in these guides from `examples/NMFS/pifsc_swordfish/`, unless a
guide explicitly specifies another location. See the [project README](../README.md)
for the build and quick start.

## Production workflow

- [Native JSON input, SS3 conversion and command-line options](NATIVE_INPUT.md)
- [Reference-point uncertainty](REFERENCE_POINT_UNCERTAINTY.md)
- [Reproducible native release bundles](NATIVE_RELEASE.md)
- [Performance and fixed-biology caching](PERFORMANCE.md)

## Model and scientific comparisons

- [Model schema, equations and seasonal dynamics](MODEL_AND_COMPARISON.md)
- [Fit diagnostics and management benchmarks](FIT_DIAGNOSTICS.md)
- [Starting from original control inputs](CONTROL_START.md)
- [Starting from an SS3 parameter snapshot](SS_PAR_START.md)
- [Dynamic unfished biomass and projections](PROJECTIONS.md)
- [Projection parameter uncertainty](PROJECTION_UNCERTAINTY.md)

## Validation and developer references

- [Assessment workflows and comparison history](ASSESSMENT_WORKFLOW.md)
- [Operational release contract](RELEASE_READINESS.md)
- [Searchable API reference](API_REFERENCE.md)
- [Independent Python conversion oracle](../tests/support/README.md)

## Code documentation and regeneration

The complete PDF is [Operations and API Guide](../output/pdf/swordfish_operations_and_api_guide.pdf).
The searchable text companion is [API_REFERENCE.md](API_REFERENCE.md); the
machine-readable inventory is `api_inventory.json`. Named declarations cover
case C++ headers/executables/tests, Python functions/classes (including nested
helpers), and R plotting functions. Anonymous lambdas, vendor internals and
generated factory bodies are explicitly excluded. The guide separately documents
the changed shared Quadra functions.

C++ library and executable definitions have Doxygen-style comments. Python
functions/classes have docstrings. Typed signatures, class fields, source
locations, module units, scientific preconditions, failure behavior, output
schemas and workflow responsibilities are covered by the guide and source.

Regenerate in this order after a source change:

```
make prepare
python3 scripts/document_api.py
/path/to/docs-environment/bin/python scripts/build_guide.py --release runs/swordfish-release-001
```

`document_api.py` requires Clang's AST JSON interface and generated input headers;
on macOS it defaults to `/usr/bin/clang++`, avoiding mismatched Homebrew SDKs.
Use `--compiler clang++` explicitly on another installation. API descriptions are
curated in `api_descriptions.json` and `python_descriptions.json`; existing Python
docstrings take precedence. Missing descriptions fail generation.

The optional `--annotate` mode inserts missing source comments/docstrings; review
the diff and run syntax/tests afterwards, then regenerate the inventory without
that flag to refresh source locations. Keep generated files separate from manual
source edits. Install PDF-only dependencies from `requirements.txt` in a virtual
environment; they are not required for fitting, release acceptance or base-R
scientific figures. Render and inspect the final PDF after regeneration.
