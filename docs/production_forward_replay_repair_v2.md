# Production forward replay repair v2

This patch matches the actual `had_quadra.hpp` structure:

```cpp
typedef unsigned int VertexId;
```

It adds graph-owned primal values and explicit value updates:

```cpp
quadra::set_value(x, 3.0);
scope.forward();
scope.zero_adjoints();
scope.backward(y);
```

Initial validated target:

```cpp
y = x * x + 1.0
```

Expected:

```text
x = 2 -> y = 5,  dy/dx = 4
x = 3 -> y = 10, dy/dx = 6
```
