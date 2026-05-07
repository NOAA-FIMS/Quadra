# Quadra

Quadra is an experimental C++ framework for sparse mixed-effects inference, Laplace approximation, and fisheries-oriented statistical modeling workflows.

It extends the original `had.h` reverse-mode automatic differentiation approach with graph-aware sparse Hessians, exact directional Hessian propagation, implicit differentiation, and scalable Laplace machinery.

## Current structure

```text
quadra/
├── quadra.hpp
├── model/
│   └── parameter.hpp
├── core/
│   ├── autodiff/
│   │   ├── autodiff.hpp
│   │   ├── adgraph.cpp
│   │   └── had/had_quadra.h
│   ├── laplace/
│   │   ├── laplace.hpp
│   │   └── evaluation.hpp
│   ├── optimizer/
│   │   └── optimizer.hpp
│   ├── inference/
│   │   └── covariance.hpp
│   ├── sparse/
│   │   └── factorization.hpp
│   └── utilities/
│       └── parallel.hpp
├── tests/
└── examples/
```

## AD graph definition

The AD graph is defined once in:

```cpp
// core/autodiff/adgraph.cpp
#include "autodiff.hpp"

DECLARE_ADGRAPH()
```

Do not put `DECLARE_ADGRAPH()` in headers.

## Build

```bash
make
```

Run tests:

```bash
make run-tests
```

Validate exact Hdot against finite differences:

```bash
make validate-hdot
```

## LaplaceOptions

```cpp
quadra::LaplaceOptions opts;
opts.use_hutchinson_trace = true;
opts.hutchinson_probes = 16;
opts.jitter_initial = 1e-12;
opts.hessian_drop_tol = 0.0;

auto fit = quadra::optimize_lbfgs(model, params, opts);
```

## Current status

Quadra currently supports:

- sparse random-effect Hessian extraction,
- graph-based Hessian sparsity discovery,
- Laplace approximation,
- implicit differentiation of optimized random effects,
- exact directional `Hdot`,
- Hutchinson trace estimation,
- adaptive sparse factorization stabilization,
- and LBFGS optimization with graceful near-convergence line-search handling.

Next milestones:

- covariance estimation,
- `sdreport`-like reporting,
- deterministic selected-inverse trace computation,
- richer fisheries likelihood examples,
- and benchmarking against TMB/FIMS workflows.
