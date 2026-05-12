# Quadra zero_adjoints support

This patch adds:

```cpp
tape.zero_adjoints();
scope.zero_adjoints();
quadra::zero_adjoints(tape);
quadra::zero_adjoints(graph);
had::ZeroAdjoints(graph);
had::ZeroAdjoints();
```

It resets accumulated reverse-mode adjoints and Hessian accumulator state without
clearing graph structure or local derivative metadata.

This fixes repeated reverse sweeps accumulating from `4` to `8`.

This is not full reusable tape support yet. Forward replay/value recomputation
is still a separate future step.
