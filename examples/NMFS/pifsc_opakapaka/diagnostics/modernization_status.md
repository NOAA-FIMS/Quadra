# opakapaka Quadra Modernization Status

This scaffold was generated as part of the Functional Analysis v1 cleanup.

## Intended layout

- `model/` — biological/model structure only
- `data/` — data row structures and loading
- `reports/` — text, CSV, and markdown report writers
- `diagnostics/` — example-specific diagnostic glue
- `quadra/` — minimal driver executable

## Next step

Move model-specific code out of the driver and wire this example to the shared
Quadra Functional Analysis report API used by the Pollock showcase.
