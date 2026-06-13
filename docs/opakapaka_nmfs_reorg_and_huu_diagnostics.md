# Opakapaka NMFS Reorganization and Huu Diagnostic Cleanup

## Status

Completed: June 2026

The PIFSC Opakapaka assessment-style example was moved under the NMFS
assessment examples directory and its final random-effect Hessian diagnostics
were corrected.

## Directory Reorganization

The Opakapaka example was moved from:

```text
examples/pifsc_opakapaka
```

to:

```text
examples/NMFS/pifsc_opakapaka
```

This keeps fisheries assessment applications separate from smaller framework
examples.

The NMFS examples directory now contains assessment-oriented examples such as:

```text
examples/NMFS/sefsc_red_snapper
examples/NMFS/pifsc_opakapaka
```

## Build Path Updates

After the move, relative include paths were updated because the example is now
one directory deeper.

For example, includes of the form:

```cpp
#include "../../../core/..."
```

were updated to:

```cpp
#include "../../../../core/..."
```

The Opakapaka executable is built from:

```bash
clang++ -std=c++17 -g -I"external/eigen/" \
  examples/NMFS/pifsc_opakapaka/quadra/opakapaka_projection.cpp \
  examples/NMFS/pifsc_opakapaka/quadra/opakapaka_adgraph_global.cpp \
  -o build/examples/pifsc_opakapaka
```

## Diagnostic Issue

After the move, the example built and ran, but the optimizer structure report
showed stale metadata:

```text
random effects     0
pattern available  no
detected structure unknown
Hessian nonzeros   0
```

This was inconsistent with the actual Laplace evaluation, which reported:

```text
Quadra: Discovering Hessian pattern from AD graph for 20 random variables ...
Quadra: Model structure aware now => Hessian pattern has 58 entries.
```

## Root Cause

The Opakapaka example can fall back to a local safeguarded one-dimensional
`log_q` polish after an L-BFGS line-search stall. That fallback returned a valid
fit and valid random effects, but it did not preserve the optimizer pattern
metadata in `fit.pattern`.

As a result, the final report was reading stale metadata even though the fitted
random-effect vector was present.

## Fix

The example now reconstructs the final random-effect Hessian after fitting:

```cpp
const Eigen::SparseMatrix<double> Huu_final =
    compute_final_random_effect_hessian(model, params, opts, fit);
```

That final Hessian is reused for:

- optimizer structure diagnostics
- Hessian nonzero reporting
- random-effect uncertainty output

This avoids relying on stale `fit.pattern` metadata when the fallback path was
used.

## Validation

After the fix, the Opakapaka example reported:

```text
random effects     20
pattern available  yes
detected structure sparse
Laplace backend    final Huu reconstruction
random solver      Laplace mode solve
Hessian nonzeros   58
```

The example also completed the fit and projection workflow and wrote outputs to:

```text
examples/NMFS/pifsc_opakapaka/outputs
```

## Remaining Note

The example still uses a local safeguarded `log_q` fallback after an L-BFGS
line-search stall:

```text
L-BFGS line-search stall detected in Opakapaka example.
Using local safeguarded one-dimensional log_q fallback.
```

This is an optimizer robustness issue, not a structural diagnostics or
uncertainty-reporting issue. The final polished fit reports a near-zero gradient
and coherent output.

Future work can replace the local fallback with a more general optimizer
robustness improvement.
