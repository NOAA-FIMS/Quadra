# Quadra core hardening patch

This patch starts splitting infrastructure out of `laplace.hpp`.

## New headers

```text
core/sparse/factorization.hpp
core/sparse/trace.hpp
core/laplace/options.hpp
```

## Equations

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

Implicit random-effect sensitivity:

```math
\frac{d\hat u}{d\theta_i}
=
-H_{uu}^{-1}H_{u\theta_i}
```

Log-determinant derivative:

```math
\frac{d}{d\theta_i}\log\det H_{uu}
=
\operatorname{tr}
\left(
H_{uu}^{-1}
\frac{dH_{uu}}{d\theta_i}
\right)
```

Hutchinson estimator:

```math
\operatorname{tr}(A) =
\mathbb{E}_z[z^\top A z]
```

## Integration note

This patch is intentionally conservative. It adds modular headers but does not require immediately deleting the existing helper implementations from `laplace.hpp`.

Next step:
- include these headers from `laplace.hpp`,
- replace duplicate helper bodies with calls into `core/sparse`,
- then move Hdot/random-effect solve code into separate Laplace modules.
