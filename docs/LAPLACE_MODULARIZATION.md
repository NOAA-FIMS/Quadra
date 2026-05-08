# Laplace modularization patch

This patch updates `core/laplace/laplace.hpp` to use the new modular sparse infrastructure:

```cpp
#include "../sparse/factorization.hpp"
#include "../sparse/trace.hpp"
#include "options.hpp"
```

## What changed

- `LaplaceOptions` now lives in `core/laplace/options.hpp`.
- adaptive sparse factorization now lives in `core/sparse/factorization.hpp`.
- trace evaluation now lives in `core/sparse/trace.hpp`.
- `laplace.hpp` keeps a small compatibility wrapper:

```cpp
logdet_directional_derivative_from_hdot(...)
```

which forwards to:

```cpp
trace_hinv_hdot(...)
```

## Key equations

Laplace approximation:

```math
\tilde f(\theta)
=
f(\theta, \hat u(\theta))
+
\frac{1}{2}\log\det H_{uu}
-
\frac{n_u}{2}\log(2\pi)
```

Log-determinant derivative:

```math
\frac{d}{d\theta_i}\log\det H_{uu}
=
\operatorname{tr}
\left(
H_{uu}^{-1}
D H_{uu}[e_i, d\hat u/d\theta_i]
\right)
```

## How to apply

Copy the patched `core/laplace/laplace.hpp` over your existing file after adding:

```text
core/sparse/factorization.hpp
core/sparse/trace.hpp
core/laplace/options.hpp
```

from the previous core-hardening patch.
