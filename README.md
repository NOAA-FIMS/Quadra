# Quadra

A high-performance C++ framework for mixed-effects modeling using Laplace approximation and automatic differentiation.

## Overview

**Quadra** is designed for fitting statistical models with both fixed and random effects, with a focus on:

- Efficient likelihood evaluation
- Scalable mixed-effects inference
- Laplace approximation for marginal likelihoods
- Gradient-based optimization via automatic differentiation

The framework is intended for applications in quantitative modeling, including fisheries stock assessment, hierarchical models, and other structured statistical systems.

## Features

- Mixed-effects modeling (fixed + random effects)
- Laplace approximation for integrating over random effects
- Automatic differentiation for exact gradients
- High-performance C++ backend
- Designed for extensibility and custom model definitions

## Build Instructions

### Requirements

- C++17-compatible compiler (e.g., clang++, g++)
- Make

### Build

```bash
make
```

### Run

```bash
make run
```

### Clean

```bash
make clean
```

## Project Structure

```
.
├── core/
│   └── eigen/        # Eigen headers (linear algebra)
├── main.cpp          # Example usage / entry point
├── Makefile          # Build configuration
└── README.md
```

## Example

A simple mixed-effects model:

```cpp
struct Model {
    double y;

    template <typename T>
    T operator()(const std::vector<T>& p) const {
        T mu = p[0];  // fixed effect
        T u  = p[1];  // random effect

        T nll = 0.0;

        // likelihood
        nll += 0.5 * (y - (mu + u)) * (y - (mu + u));

        // random effect prior
        nll += 0.5 * u * u;

        return nll;
    }
};
```

## Design Philosophy

- Explicit over implicit: models are defined directly in C++
- Performance-first: minimal abstraction overhead
- Mathematically grounded: aligns closely with statistical formulation
- Composable: supports building more complex hierarchical systems

## Roadmap

- Sparse precision matrix support
- Advanced optimization routines (e.g., L-BFGS, trust region)
- Extended distribution library
- R/Python interface bindings
- Parallel evaluation

## Related Work

This framework is conceptually similar to:

- TMB (Template Model Builder)
- Stan

…but aims to provide a lightweight, transparent alternative with fine-grained control over model structure and computation.

## License

Public Domain

## Author

Matthew Supernaw
