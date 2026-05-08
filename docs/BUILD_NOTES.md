# Build notes

The AD graph should still be defined once:

```cpp
// core/autodiff/adgraph.cpp
#include "autodiff.hpp"

DECLARE_ADGRAPH()
```

Link `core/autodiff/adgraph.cpp` into tests and examples:

```make
CORE_SRC = core/autodiff/adgraph.cpp
```

Use project-root relative includes such as:

```cpp
#include "quadra.hpp"
```

or subsystem headers:

```cpp
#include "core/inference/inference.hpp"
```
