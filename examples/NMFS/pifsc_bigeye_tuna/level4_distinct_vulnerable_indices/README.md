# Level 4: Distinct Vulnerable-Biomass Indices

Level 3 added fleet-specific selectivity and fleet-vulnerable-biomass index
predictions. The objective improved, but fixed-effect geometry still showed
weak abundance/F/q directions.

The index ratio inspection suggested the synthetic fleet indices were still too
similar. They carried almost the same abundance trend, so they behaved like
fleet labels attached to a pooled index.

## Model change

Compared with Level 3:

- keep the same model structure
- regenerate synthetic fleet indices from deliberately different vulnerable
  biomass signals
- longline index emphasizes older fish
- purse seine index emphasizes younger fish

## Scientific hypothesis

If the main problem was that fleet indices were effectively pooled, then
distinct vulnerable-biomass index trends should improve separability among:

- R0
- Fbar
- q_longline
- q_purse_seine
- fleet-specific selectivity

## Decision rule

If distinct fleet index trajectories improve geometry, the workflow has
identified a data-information bottleneck rather than an optimizer failure.
