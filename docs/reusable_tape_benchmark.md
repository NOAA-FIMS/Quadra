# Reusable tape benchmark

This benchmark compares repeated gradient evaluations for a random-intercept
joint objective using two modes:

1. Rebuild the AD tape for every evaluation.
2. Build the AD tape once, then use:

```cpp
quadra::set_value(mu, new_mu);
quadra::set_value(u, new_u);
scope.forward();
scope.zero_adjoints();
scope.backward(nll);
```

Run:

```bash
make tests/test_reusable_tape_random_intercept_gradient
./tests/test_reusable_tape_random_intercept_gradient

make benchmark-reusable-tape
```
