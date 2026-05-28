# Quadra Exact Laplace Gradient Architecture

## Overview

Quadra is a modern mixed-effects inference framework focused on exact directional differentiation, structure-aware sparse computations, and reusable automatic differentiation infrastructure.

The framework originated from experimentation with reverse-mode automatic differentiation and Laplace approximation methods, but has evolved toward a broader goal:

> Build a modular inference engine where model structure, derivative structure, and computational reuse are first-class concepts.

Rather than treating the Laplace approximation as a black-box numerical procedure, Quadra exposes and optimizes the mathematical structure underlying mixed-effects inference.

---

# Core Design Philosophy

Quadra is intentionally not a clone of existing AD-based mixed-effects systems.

The goal is not simply:

```text
objective function + automatic differentiation
```

Instead, the framework is built around:

```text
structure-aware differentiation
```

This means:

- exploiting sparsity patterns,
- reusing derivative state,
- differentiating Hessian structure directly,
- minimizing repeated graph traversals,
- and separating model semantics from differentiation mechanics.

The framework emphasizes modularity and composability over monolithic modeling interfaces.

---

# Exact Laplace Gradient Stack

The Laplace approximation objective is:

```text
L(theta)
=
f(theta, u*)
+
0.5 log det(H_uu(theta, u*))
```

where:

- `theta` are fixed effects,
- `u*` are optimized random effects,
- `H_uu` is the random-effects Hessian.

Quadra computes exact directional derivatives of the Laplace objective through explicit Hessian derivative propagation rather than relying solely on finite differences.

The current architecture contains several interacting layers.

---

# Directional Hdot Propagation

A central component of the framework is exact directional Hessian differentiation:

```text
Hdot = D H_uu(theta, u*) [direction]
```

This is implemented through directional edge-pushing extensions to the `had_quadra.hpp` automatic differentiation system.

Key additions include:

- directional tangent propagation,
- directional adjoint propagation,
- directional second-order edge accumulation,
- sparse directional Hessian edge storage,
- exact directional Hessian extraction.

This avoids repeated finite-difference Hessian rebuilds and enables exact trace derivatives.

---

# Replay-Reuse AD Execution

Repeated retaping and replaying of large AD graphs can dominate runtime.

Quadra introduces replay-reuse infrastructure that separates:

- graph construction,
- directional seeding,
- Hessian extraction,
- and trace contraction.

This allows:

- reuse of graph topology,
- reuse of sparsity patterns,
- reuse of replay state,
- and reuse of factorized Hessians.

Replay-reuse becomes increasingly important for large mixed-effects systems and repeated optimizer evaluations.

---

# Active-Direction Discovery

Many fixed-effect directions produce zero or structurally sparse Hessian derivatives.

Quadra includes active-direction discovery mechanisms that identify directions that actually contribute to:

```text
D H_uu(theta, u*)
```

This allows the framework to skip unnecessary directional propagation and trace contraction work.

The resulting system becomes increasingly structure-aware as model dimensionality grows.

---

# Lazy Implicit Random-Effect Sensitivities

Traditional implementations frequently materialize the full dense matrix:

```text
du*/dtheta
```

Quadra instead supports lazy directional access:

```text
du*/dtheta_j
```

computed only when needed.

This enables:

- sparse active-direction workflows,
- reduced memory traffic,
- factorization reuse,
- and optimizer-aware directional scheduling.

The current implementation supports reusable solve providers:

```text
du*/dtheta_j = - H_uu^{-1} f_{u theta_j}
```

evaluated lazily through reusable Hessian factorizations.

---

# Cached Trace Contraction

The Laplace gradient contains trace terms of the form:

```text
trace(H^{-1} Hdot)
```

Quadra includes reusable contraction infrastructure that:

- caches sparse solves,
- reuses selected Hessian columns,
- avoids repeated dense contractions,
- and separates trace contraction from directional propagation.

This becomes increasingly important as model dimension and directional count increase.

---

# Optimizer-Aware Architecture

The framework is moving toward optimizer-aware execution.

The long-term goal is to allow optimizers to interact directly with reusable derivative infrastructure instead of repeatedly rebuilding derivative state.

This includes:

- reusable directional providers,
- reusable Hessian factorizations,
- reusable sparse solve caches,
- active-direction scheduling,
- and eventually batched directional propagation.

The architecture is designed so that future optimizers can request only the derivative information required at a particular iteration.

---

# Why Quadra Is Different

Many mixed-effects systems treat automatic differentiation as a general-purpose black box.

Quadra instead treats derivative structure itself as an optimization target.

The framework is built around:

- exact directional derivative propagation,
- structure-aware sparsity,
- replay reuse,
- implicit sensitivity methods,
- and modular derivative engines.

The resulting system is intended to scale not only numerically, but architecturally.

---

# Current Research Directions

Ongoing areas of experimentation include:

- exact third-order directional derivatives,
- sparse Hessian derivative propagation,
- batched directional edge-pushing,
- implicit trace differentiation,
- Hessian-vector and Hessian-Hessian products,
- adaptive active-direction discovery,
- and optimizer-integrated derivative scheduling.

These directions move Quadra beyond a traditional AD wrapper toward a reusable computational inference framework.

---

# Long-Term Vision

The long-term vision for Quadra is:

```text
A modular, structure-aware inference engine
for modern mixed-effects statistical computing.
```

The framework aims to provide:

- scalable derivative infrastructure,
- composable model abstractions,
- reusable computational graphs,
- exact directional differentiation,
- and modern sparse numerical workflows.

Quadra is intended to support research-grade experimentation while remaining suitable for production scientific inference systems.
