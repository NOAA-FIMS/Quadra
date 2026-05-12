# had_quadra adjoint reset probe

This diagnostic checks whether repeated calls to `scope.backward(y)` accumulate
adjoints and whether the current API exposes reset hooks such as:

- `zero_adjoints()`
- `reset_adjoints()`
- `clear_adjoints()`
- derivative reset variants

The previous mutation test produced:

```text
initial grad = 4
mutated grad = 8
```

This probe tests whether the `8` is simply the result of running backward twice
without clearing adjoints.
