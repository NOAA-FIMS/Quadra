# Cached Laplace objective

This patch adds:

```cpp
quadra::CachedLaplaceObjectiveState
quadra::evaluate_laplace_objective_cached(...)
```

The cached evaluator reuses symbolic sparse factorization for the random-effect
Hessian logdet calculation.

This is additive and does not replace the existing `evaluate_laplace_objective`.

Run:

```bash
make tests/test_laplace_objective_cached
./tests/test_laplace_objective_cached
make benchmark-cached-laplace-objective
```
