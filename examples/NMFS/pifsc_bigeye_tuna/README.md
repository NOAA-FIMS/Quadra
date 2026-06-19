# PIFSC Bigeye Tuna Diagnostic-Guided Modeling Prototype

This example is intentionally **not** a MULTIFAN-CL port.

The purpose is to build a sequence of increasingly complex highly migratory
species models and use Quadra diagnostics to decide which model features are
supported by the data.

The core workflow is:

1. Start with the simplest plausible model.
2. Fit the model.
3. Inspect convergence, curvature, identifiability, and data influence.
4. Add one biological or observation-process feature.
5. Refit and compare diagnostics.
6. Retain complexity only when the diagnostics support it.

## Planned modeling ladder

| Level | Model feature | Diagnostic question |
| --- | --- | --- |
| 0 | Single-region age-structured model | Are abundance, recruitment, and fishing mortality separable? |
| 1 | Multiple fleets | Which fleets identify abundance, selectivity, and catchability? |
| 2 | Fleet-specific composition data | Which data components identify selectivity? |
| 3 | Multiple regions | Are regional abundance trends separately identified? |
| 4 | Movement | Is movement identifiable or confounded with recruitment/q? |
| 5 | Tagging | Does tagging identify movement, reporting, or both? |

Each level should produce:

- fit summary
- objective components
- functional analysis report
- Laplace structure report
- reference points, when appropriate
- model decision notes

This prototype is synthetic and public-data-safe. It is not an official
assessment.
