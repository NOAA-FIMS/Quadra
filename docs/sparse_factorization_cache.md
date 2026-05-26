# Sparse factorization cache

This patch adds:

```cpp
quadra::SparseFactorizationCache
```

The goal is to reuse symbolic sparse factorization:

```cpp
cache.analyze_pattern(H_pattern); // once
cache.factorize(H_values);        // many times
cache.logdet();
cache.solve(b);
```

This avoids repeated:

```cpp
solver.compute(H);
```

when the Hessian sparsity pattern is stable but numeric values change.

Initial benchmark:

```bash
make benchmark-sparse-factorization-cache
```

This benchmark uses repeated tridiagonal SPD matrices, representative of AR1
random-effect Hessians.
