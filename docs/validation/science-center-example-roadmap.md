# Science Center Example Validation Roadmap

This document tracks a proposed validation suite with one representative assessment-style example from each NOAA Fisheries Science Center.

The goal is to build examples that are:
- public-data-safe or synthetic,
- reproducible,
- paired with TMB reference implementations where practical,
- documented with expected outputs,
- capable of reporting uncertainty, derived quantities, and projections.

## Proposed example set

| Science Center | Example | Status | Main validation target |
|---|---|---:|---|
| PIFSC | Opakapaka projection example | In progress | Projection validation and Level-1 uncertainty reporting |
| SEFSC | Red-snapper-style age-structured model | Scaffolded | Age structure, selectivity, recruitment deviations, projections |
| NEFSC | Groundfish/index-heavy assessment | Planned | Multiple indices, survey likelihoods, retrospective-style diagnostics |
| NWFSC | West Coast age-structured model | Planned | Age composition, selectivity, biological reference points |
| AFSC | Pollock/sablefish-style model | Planned | Recruitment deviations, state-space/random-effect scalability |
| SWFSC | CPS/tuna-style model | Planned | Time-varying dynamics, index scaling, projection scenarios |

## Shared validation requirements

Each example should eventually include:

1. Quadra implementation
2. TMB comparison implementation
3. synthetic or public-data-safe input data
4. reproducible runner
5. fit diagnostics
6. standard errors and confidence intervals
7. random-effect conditional uncertainty
8. derived quantity uncertainty
9. projection envelopes
10. comparison summary against TMB

## Recommended directory layout

```text
examples/<example_name>/
  README.md
  data/
  quadra/
  tmb/
  outputs/
  validation/
```

## Development order

1. Finish Opakapaka Level-1 uncertainty reporting.
2. Scaffold SEFSC red-snapper-style age-structured model.
3. Add minimal Quadra implementation.
4. Add TMB reference implementation.
5. Add validation summary and uncertainty outputs.
6. Repeat for the remaining science centers.
