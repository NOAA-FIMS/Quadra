# Swordfish code documentation

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
