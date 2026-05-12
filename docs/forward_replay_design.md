# Quadra forward replay design

## Motivation

Recent diagnostics established:

```text
zero_adjoints = solved
forward replay = missing
```

`x.val` can be mutated, and reverse adjoints can now be reset, but downstream
primal values remain stale because the current graph does not replay operations.

To support reusable tapes, Quadra needs the ability to:

```text
build graph once
set independent values
forward replay values
zero adjoints
reverse sweep
```

## Target API

A future high-level API could look like:

```cpp
quadra::ReusableTape tape(model, initial_params);

tape.set_values(new_params);
double objective = tape.forward();
Eigen::VectorXd gradient = tape.reverse();
```

Lower-level AD API:

```cpp
quadra::set_independent_values(graph, inputs, values);
quadra::forward(graph);
quadra::zero_adjoints(graph);
scope.backward(objective);
```

## Minimum graph metadata

Each replayable graph vertex needs enough information to recompute its primal
value from parent values.

A minimal operation code:

```cpp
enum class OpCode {
    Constant,
    Independent,
    Add,
    Subtract,
    Multiply,
    Divide,
    Exp,
    Log,
    Sqrt,
    Negate
};
```

Each vertex needs:

```cpp
OpCode op;
VertexId left;
VertexId right;
double constant;
bool has_left;
bool has_right;
```

For unary operations, only `left` is used.
For constants and independents, no parents are used.

## Forward replay rule

The graph is already constructed in topological order, so replay can be:

```cpp
for vertex in vertices:
    switch(vertex.op):
        case Independent:
            keep assigned value
        case Add:
            value = value(left) + value(right)
        case Multiply:
            value = value(left) * value(right)
        ...
```

## Why prototype first

The production `had_quadra.hpp` currently stores local derivative information
for reverse propagation, but it does not yet expose operation-level replay
metadata.

Before modifying all operator overloads, this patch adds an isolated prototype
test showing the target replay behavior for:

```text
y = x*x + 1
```

This validates the architecture independently from the existing AD graph.
