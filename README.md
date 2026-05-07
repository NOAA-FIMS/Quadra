# Quadra

Quadra is a modern C++ framework for sparse mixed-effects inference, Laplace approximation, and scalable statistical modeling.

It is designed as a clean, modular alternative to large monolithic inference systems, with an emphasis on:

- explicit APIs,
- sparse computation,
- scalable random effects,
- graph-aware automatic differentiation,
- exact directional Hessian propagation,
- and composable inference/reporting layers.

Quadra was originally motivated by fisheries stock assessment workflows, but the architecture is general and suitable for many latent-variable and hierarchical modeling problems.

---

# Core Design Philosophy

Quadra intentionally avoids:

- macro-heavy APIs,
- hidden global registries,
- implicit reporting side effects,
- opaque AD behavior,
- and tightly coupled monolithic workflows.

Instead, Quadra emphasizes:

```text
models are ordinary C++ classes
inference is explicit
reporting is explicit
uncertainty is composable
```

The framework separates:

- objective evaluation,
- optimization,
- Laplace approximation,
- covariance estimation,
- and reporting

into clean independent layers.

---

# Current Features

## Automatic Differentiation

Quadra extends the `had.h` reverse-mode automatic differentiation framework with:

- graph-aware sparse Hessian discovery,
- exact directional Hessian propagation (`Hdot`),
- sparse Hessian extraction,
- implicit differentiation,
- and scalable sparse mixed-effects workflows.

---

## Sparse Laplace Approximation

Implemented:

- sparse random-effect Hessians,
- graph-based sparsity discovery,
- sparse LDLT factorization,
- adaptive Hessian stabilization,
- implicit differentiation of optimized random effects,
- Hutchinson trace approximation,
- exact directional curvature propagation,
- and scalable random-effect optimization.

---

## Optimization

Current optimization support:

- L-BFGS optimization,
- Newton optimization for random effects,
- graceful line-search handling near convergence,
- sparse Hessian-aware workflows,
- and adaptive jitter stabilization.

---

## Covariance Estimation

Quadra includes first-pass covariance estimation for fixed effects:

```text
Cov(theta_hat) ≈ H_laplace(theta_hat)^-1
```

Current implementation:
- finite-difference Hessian of the optimized Laplace gradient,
- adaptive SPD stabilization,
- standard error storage directly on parameters.

Future versions will support:
- exact Laplace Hessians,
- selected sparse inverse methods,
- derived quantity covariance,
- random-effect covariance,
- and sdreport-style workflows.

---

# Report Registry System

Quadra uses an explicit report registry rather than macro-based reporting.

## Example

```cpp
template <typename T, typename ReportLike>
void report(const std::vector<T>& p, ReportLike& out) const {

    T ssb = ...;
    T msy = ...;

    out.estimate("biomass/SSB", ssb);

    out.estimate(
        "reference_points/MSY",
        msy);

    out.add(
        "diagnostics/max_gradient",
        max_gradient);
}
```

---

## Why Explicit Reporting?

Optimization never touches reporting.

Reporting is evaluated only post-fit:

```cpp
auto fit =
    optimize_lbfgs(model, params, laplace_options);

auto covariance =
    estimate_fixed_covariance(
        model,
        params,
        fit);

auto report =
    evaluate_report_with_uncertainty(
        model,
        params,
        covariance);
```

This avoids:
- report registration during optimization,
- hidden side effects,
- and macro-driven AD behavior.

---

## Hierarchical Report Paths

Reports support hierarchical names:

```cpp
out.estimate("biomass/SSB_2025", ssb);
out.estimate("reference_points/MSY", msy);
out.add("diagnostics/max_gradient", max_grad);
```

Groups are inferred automatically from paths.

---

## Report Metadata

Reports support metadata:

```cpp
quadra::ReportMetadata meta;

meta.units = "metric tons";
meta.description = "Spawning stock biomass";

out.estimate(
    "biomass/SSB",
    ssb,
    meta);
```

---

## CSV Export

```cpp
report.to_csv("report.csv");
```

CSV output includes:

```text
path
group
name
value
std_error
estimate_uncertainty
units
description
```

---

# Project Structure

```text
quadra/
├── quadra.hpp
├── model/
│   └── parameter.hpp
├── core/
│   ├── autodiff/
│   │   ├── autodiff.hpp
│   │   ├── adgraph.cpp
│   │   └── had/
│   │       └── had_quadra.h
│   │
│   ├── optimizer/
│   │   └── optimizer.hpp
│   │
│   ├── laplace/
│   │   ├── laplace.hpp
│   │   └── evaluation.hpp
│   │
│   ├── inference/
│   │   ├── covariance.hpp
│   │   └── report.hpp
│   │
│   ├── sparse/
│   │   └── factorization.hpp
│   │
│   └── utilities/
│       └── parallel.hpp
│
├── tests/
├── examples/
└── docs/
```

---

# AD Graph Definition

The AD graph should be defined exactly once:

```cpp
// core/autodiff/adgraph.cpp

#include "autodiff.hpp"

DECLARE_ADGRAPH()
```

Do not place `DECLARE_ADGRAPH()` in headers.

---

# Build

## Build all tests and examples

```bash
make
```

## Run tests

```bash
make run-tests
```

## Validate exact Hdot propagation

```bash
make validate-hdot
```

---

# Example Model Structure

```cpp
struct MyModel {

    template <typename T>
    T objective(const std::vector<T>& p) const {

        T nll = 0.0;

        return nll;
    }

    template <typename T>
    T operator()(const std::vector<T>& p) const {
        return objective(p);
    }

    template <typename T, typename ReportLike>
    void report(
        const std::vector<T>& p,
        ReportLike& out) const {

        T ssb = ...;

        out.estimate("biomass/SSB", ssb);
    }
};
```

---

# Current Development Priorities

## Near-Term

- exact Laplace Hessians,
- selected sparse inverse methods,
- faster report uncertainty propagation,
- random-effect covariance estimation,
- derived quantity covariance,
- fisheries-specific examples,
- AR1 and state-space workflows,
- and benchmarking against TMB/FIMS workflows.

## Longer-Term

- sparse sdreport equivalent,
- simulation utilities,
- MCMC workflows,
- JSON export,
- parallel sparse factorization,
- and distributed inference workflows.

---

# Acknowledgements

Quadra builds upon ideas and techniques developed across the scientific computing, automatic differentiation, and mixed-effects modeling communities.

Special acknowledgement goes to:

- Bachi Li for the original `had.h` reverse-mode automatic differentiation framework and higher-order derivative work that inspired Quadra's AD foundation.
- The developers of Eigen for providing the linear algebra infrastructure used throughout the project.
- The developers of TMB (Template Model Builder), whose innovations in sparse Laplace approximation and mixed-effects inference strongly influenced the broader design space explored by Quadra.
- NOAA Fisheries and the Fisheries Integrated Modeling System (FIMS) project for motivating scalable and maintainable mixed-effects workflows in fisheries stock assessment.

Quadra is ultimately an experiment in building a cleaner, more modular inference architecture for large scientific modeling systems.

---

# Status

Quadra is currently an experimental framework under active development.

The architecture is already capable of:
- scalable sparse mixed-effects optimization,
- Laplace approximation,
- covariance estimation,
- and uncertainty propagation for derived quantities.

The long-term goal is a modern, modular inference framework suitable for large scientific modeling systems without the dependency and architectural complexity often associated with legacy mixed-effects frameworks.
