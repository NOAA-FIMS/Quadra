# Quadra

A lightweight mixed-effects inference framework for modern scientific computing.

Quadra is an experimental C++ framework for sparse mixed-effects optimization and Laplace approximation, designed as a clean, modular alternative to heavier AD-based systems.

The project began as an extension of Bachi Li's `had.h` reverse-mode automatic differentiation library and has since evolved into a broader framework for:

* sparse mixed-effects inference,
* Laplace approximation,
* exact directional Hessian propagation,
* graph-aware Hessian sparsity discovery,
* implicit differentiation,
* scalable random-effect optimization,
* and fisheries-oriented statistical modeling workflows.

Quadra is being developed with long-term goals of supporting large-scale stock assessment workflows and eventually providing a lightweight alternative to Template Model Builder (TMB) for systems such as FIMS.

---

# Current Features

## Automatic Differentiation

Quadra extends the original `had.h` reverse-mode AD framework with additional capabilities used for mixed-effects inference:

* reverse-mode automatic differentiation,
* Hessian extraction,
* directional Hessian propagation,
* sparse graph-aware derivative workflows,
* and experimental third-order derivative infrastructure.

## Sparse Mixed-Effects Inference

Current mixed-effects features include:

* Laplace approximation,
* sparse random-effect Hessian extraction,
* graph-based Hessian sparsity discovery,
* exact directional Hessian propagation (`Hdot`),
* implicit differentiation for `du*/dtheta`,
* adaptive sparse factorization stabilization,
* Hutchinson trace estimation,
* sparse LDLT factorization,
* and scalable random-effect optimization.

## Optimization

Quadra currently uses:

* LBFGS++ for outer optimization,
* sparse Newton solves for random effects,
* adaptive line-search recovery near convergence,
* and sparse factorization reuse.

## Fisheries-Oriented Modeling

The framework is being developed with fisheries stock assessment workflows in mind.

Current examples include:

* CPUE/index models,
* latent abundance process models,
* random year effects,
* selectivity random walks,
* Gaussian random walks,
* Poisson random effects,
* and curvature-dependent validation models.

---

# Project Structure

```text
quadra/
├── core/
│   ├── autodiff/
│   ├── laplace/
│   ├── optimizer/
│   ├── sparse/
│   └── utilities/
│
├── tests/
│   ├── test_curvature_depends_on_theta.cpp
│   ├── test_poisson_random_effect.cpp
│   ├── test_ar1_random_walk.cpp
│   └── test_hdot_validation.cpp
│
├── examples/
│   ├── fisheries_random_year_effects.cpp
│   ├── fisheries_age_selectivity_random_walk.cpp
│   └── fisheries_index_cpue_laplace.cpp
│
├── docs/
└── README.md
```

---

# Build Requirements

## Dependencies

Quadra currently depends on:

* C++17 or newer,
* Eigen,
* LBFGS++,
* and the modified `had.h` implementation.

## Example Build

```bash
clang++ -std=c++17 -O3 -flto \
  -Icore/eigen \
  -o examples/fisheries_random_year_effects \
  examples/fisheries_random_year_effects.cpp \
  core/adgraph.cpp
```

---

# AD Graph Definition

Quadra uses a single global/thread-local AD graph definition.

Create:

```cpp
// core/adgraph.cpp
#include "autodiff.hpp"

DECLARE_ADGRAPH()
```

This source file should be linked into all tests/examples/executables.

Avoid placing `DECLARE_ADGRAPH()` directly inside headers.

---

# LaplaceOptions

Quadra exposes configurable Laplace approximation settings through:

```cpp
quadra::LaplaceOptions
```

Example:

```cpp
quadra::LaplaceOptions opts;

opts.use_hutchinson_trace = true;
opts.hutchinson_probes = 16;
opts.jitter_initial = 1e-12;
opts.jitter_max_attempts = 12;
opts.hessian_drop_tol = 0.0;

auto fit = quadra::optimize_lbfgs(model, params, opts);
```

## Available Options

| Option                 | Description                                |
| ---------------------- | ------------------------------------------ |
| `use_hutchinson_trace` | Use stochastic Hutchinson trace estimation |
| `hutchinson_probes`    | Number of Hutchinson probes                |
| `hutchinson_seed`      | RNG seed for stochastic trace estimation   |
| `jitter_initial`       | Initial diagonal stabilization jitter      |
| `jitter_max_attempts`  | Maximum adaptive jitter attempts           |
| `validate_hdot`        | Enable exact-vs-FD Hdot validation         |
| `hessian_drop_tol`     | Sparse Hessian entry dropping threshold    |

---

# Exact Directional Hessian Propagation

One of the major architectural features of Quadra is exact directional Hessian propagation:

```text
Hdot = D H_uu [e_i, du*/dtheta_i]
```

This replaces older finite-difference Hessian rebuild workflows and significantly improves scalability for large random-effect systems.

The current implementation supports:

* exact directional derivative propagation,
* sparse Hessian extraction,
* and stochastic or deterministic trace evaluation.

---

# Testing

## Build All Tests

```bash
make
```

## Run Tests

```bash
make run-tests
```

## Validate Exact Hdot

```bash
make validate-hdot
```

This enables:

```bash
-DQUADRA_VALIDATE_HDOT
```

which compares exact directional Hessian propagation against the finite-difference validation path on small systems.

---

# Current Development Goals

Short-term development goals include:

* exact deterministic sparse trace evaluation,
* selected inverse methods,
* sparse Cholesky reuse,
* covariance estimation (`sdreport`-like functionality),
* report quantities,
* parallel Hessian extraction,
* sparse higher-order derivatives,
* and fisheries-specific modeling workflows.

Long-term goals include:

* scalable stock assessment inference,
* integration with fisheries modeling systems,
* improved sparse mixed-effects infrastructure,
* and a lightweight alternative to TMB-style workflows.

---

# Status

Quadra is currently experimental research software.

The framework is under active development and APIs may change frequently.

---

# Acknowledgements

Quadra builds upon and extends concepts from:

* Bachi Li's `had.h` automatic differentiation framework,
* sparse mixed-effects inference methods,
* Laplace approximation techniques,
* and scientific computing workflows used in fisheries stock assessment systems.

The current implementation includes substantial extensions and modifications developed specifically to support the Quadra framework and sparse mixed-effects inference workflows.
