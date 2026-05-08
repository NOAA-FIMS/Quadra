# Quadra modularization plan

This patch establishes stable module boundaries without moving all implementation code yet.

## Public includes

```cpp
#include "quadra.hpp"
```

## Internal umbrella headers

```text
core/sparse/sparse.hpp
core/inference/inference.hpp
core/laplace/laplace_all.hpp
```

## Laplace module boundaries

```text
core/laplace/options.hpp
core/laplace/evaluation.hpp
core/laplace/random_effects.hpp
core/laplace/implicit_diff.hpp
core/laplace/hdot.hpp
core/laplace/laplace.hpp
```

## Next refactor sequence

1. Move `solve_random_effects_laplace(...)` into `random_effects.hpp`.
2. Move `implicit_du_dtheta_*` into `implicit_diff.hpp`.
3. Move `random_hessian_directional_exact(...)` into `hdot.hpp`.
4. Keep `laplace.hpp` as the high-level orchestrator.
5. Add unit tests after each move.

## Architectural rule

A module should know only what it must know.

For example:
- `core/sparse` should not know about models.
- `core/inference` should not know fisheries logic.
- `core/laplace` should not own reporting.
- `model` should not own optimizer behavior.
