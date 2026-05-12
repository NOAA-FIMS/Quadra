# Reusable tape plan for Quadra / had_quadra

## Goal

Avoid repeatedly rebuilding the same AD graph during repeated Laplace evaluations.

Current pattern:

```text
build tape -> evaluate model -> reverse sweep -> discard tape
```

Target pattern:

```text
build tape once -> update independent values -> forward/value pass -> reverse sweep
```

## Required capability

The Quadra AD layer needs an explicit way to update the values of independent
variables attached to an existing graph.

A desirable API might look like:

```cpp
quadra::ReusableTape tape(model, initial_params);

tape.set_values(new_params);
double objective = tape.forward();
Eigen::VectorXd gradient = tape.reverse();
```

or at the lower level:

```cpp
auto inputs = quadra::make_independent(tape, values);
quadra::set_independent_values(inputs, new_values);
quadra::forward(tape.graph);
quadra::backward(objective);
```

## Why this matters

Recent benchmarks indicate the dominant cost is not raw double evaluation or
tape construction alone, but repeated AD reverse/Hessian work across repeated
Laplace and LBFGS evaluations. Reusable tapes are the next likely optimization
path.

## Current patch

This patch adds a capability probe:

```bash
make tests/test_had_quadra_reusable_tape_probe
./tests/test_had_quadra_reusable_tape_probe
```

The result tells us whether mutable independent-value support is already exposed
or whether the next patch should add it to `had_quadra.hpp` / `autodiff.hpp`.
