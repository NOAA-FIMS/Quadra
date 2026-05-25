# Quadra

Quadra is an experimental mixed-effects inference framework focused on:

- reverse-mode automatic differentiation
- Laplace approximation
- sparse Hessian methods
- implicit differentiation
- scalable scientific inference workflows

The project is oriented toward computational statistics, state-space models,
and large-scale scientific modeling applications such as fisheries stock
assessment.

Current development emphasizes:

- transparent numerical workflows
- reproducible benchmarking
- sparse linear algebra performance
- modular inference architecture
- computational profiling and scaling analysis

---

# Core Capabilities

Current implemented capabilities include:

- reverse-mode automatic differentiation
- exact and approximate gradients
- Laplace approximation
- random effects optimization
- sparse Hessian construction
- sparse factorization reuse
- implicit derivative calculations
- profiled derived quantity uncertainty
- delta-method utilities
- benchmark and scaling infrastructure
- benchmark normalization and plotting
- CI benchmark workflows

---

# Repository Structure

```text
core/
    autodiff/
    laplace/
    inference/
    model/
    optimizer/

examples/
    simple/
    big/

benchmarks/
    analysis/
    comparisons/
    outputs/
    normalized/

tests/
```

---

# Building

## Requirements

Recommended:

- C++17 compiler
- Eigen
- LBFGSpp
- GNU Make

Optional:

- R
- TMB

## Build benchmarks and examples

```bash
make
```

---

# Running Tests

Run all tests:

```bash
make test
```

Run selected tests:

```bash
make test-laplace-implicit-workspace
make test-laplace-profiled-derived-report
make test-random-intercept
```

---

# Running Examples

## Simple random intercept example

```bash
./examples/simple/random_intercept_model
```

## Catch-at-age Laplace example

```bash
./examples/big/catch_at_age_laplace
```

---

# Benchmarking

## Random intercept scaling

```bash
make benchmark-random-intercept
```

## State-space scaling

```bash
make benchmark-state-space
```

## Quadra vs TMB comparison scaffold

Quadra benchmark:

```bash
make benchmark-quadra-tmb-random-intercept-quadra
```

TMB benchmark:

```bash
make benchmark-quadra-tmb-random-intercept-tmb
```

---

# Benchmark Analysis

Normalize benchmark outputs:

```bash
make benchmark-normalize-all
```

Generate scaling plots:

```bash
make benchmark-plot-random-intercept
```

---

# Continuous Integration

GitHub Actions workflows currently provide:

- contract compilation checks
- benchmark execution
- RSS logging
- benchmark normalization
- benchmark artifact uploads
- scaling plot generation

---

# Current V1 Focus

The current V1 effort is focused on:

- stable mixed-effects workflows
- reusable sparse factorization infrastructure
- implicit differentiation utilities
- scalable state-space inference
- reproducible benchmark infrastructure
- comparative inference benchmarking

---

# Project Status

Quadra is currently experimental and under active development.

Interfaces, APIs, benchmark structure, and numerical methods may evolve rapidly
during the V1 development cycle.
