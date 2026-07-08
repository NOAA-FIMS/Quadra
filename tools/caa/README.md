# CAA Toolchain

The CAA Toolchain transforms declarative assessment descriptions into generated documentation, validation reports, and eventually runtime components.

The runtime does not depend on the toolchain.

## Inputs

- `caa_manifest.yml`
- `package.meta`

## Outputs

- package catalog
- package registry
- architecture summaries
- validation reports
- generated assessment-cycle code

## Principles

- Metadata is the single source of truth.
- Generated files are never edited by hand.
- Scientific code is independent of generators.
- Every generated artifact should be reproducible.

## Flow

```text
Manifest + Metadata
        |
        v
   CAA Toolchain
        |
        v
Generated Artifacts
        |
        v
 AssessmentCycle
        |
        v
 Quadra Solver
```
