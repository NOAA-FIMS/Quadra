# Quadra validation suite plan

## Core validation tests

```text
test_hdot_exact_vs_fd.cpp
test_implicit_diff.cpp
test_laplace_gradient_consistency.cpp
test_covariance.cpp
test_report_registry.cpp
test_sparse_factorization.cpp
test_trace_estimators.cpp
```

## Scaling tests

```text
test_large_random_effects_1k.cpp
test_large_random_effects_10k.cpp
test_large_random_effects_100k.cpp
```

These should be opt-in rather than part of the default test suite.

## Mathematical checks

### Hdot

```math
\dot H
=
D H_{uu}[v]
```

Compare exact directional propagation against central finite difference.

### Implicit differentiation

```math
\frac{d\hat u}{d\theta_i}
=
-H_{uu}^{-1}H_{u\theta_i}
```

Compare against finite-difference changes in optimized random effects.

### Laplace gradient

```math
\nabla_\theta \tilde f(\theta)
```

Compare analytic implementation against central finite differences of the Laplace objective on small models.

### Covariance

```math
\operatorname{Cov}(\hat\theta)
\approx
H_{\tilde f}^{-1}
```

Validate against analytic toy models where possible.
