# Forward replay math ops

This patch extends production replay support to:

- division
- division by constants
- constants divided by AD variables
- `exp`
- `log`
- `sqrt`

Validation expression:

```cpp
y = exp(x) + log(x) + sqrt(x) + 10.0 / x + x / 2.0;
```

Replay flow:

```cpp
quadra::set_value(x, 3.0);
scope.forward();
scope.zero_adjoints();
scope.backward(y);
```
