# Quadra V1 Workflow

## What Quadra Is

Quadra is a modern mixed-effects inference framework focused on:

- Laplace approximation
- automatic differentiation
- sparse Hessian methods
- implicit differentiation
- derived quantity uncertainty
- scalable scientific modeling workflows

The framework is designed around composable inference primitives that can
support both small statistical models and large production-scale systems such
as catch-at-age stock assessment models.

---

# Core Concepts

## Fixed Effects

Fixed effects are the primary parameters estimated directly by optimization.

Examples:

- growth parameters
- observation variance
- selectivity parameters
- regression coefficients

In Quadra, fixed effects are represented by:

```cpp
false // random = false
```

during parameter registration.

---

## Random Effects

Random effects are latent variables integrated using the Laplace approximation.

Examples:

- annual deviations
- latent trajectories
- state-space processes
- random intercepts

In Quadra, random effects are represented by:

```cpp
true // random = true
```

during parameter registration.

---

# The Quadra Inference Workflow

The standard V1 inference workflow is:

```text
Define model
↓
Register parameters
↓
Partition fixed/random effects
↓
Optimize Laplace objective
↓
Build Laplace implicit workspace
↓
Estimate fixed-effect covariance
↓
Compute profiled derived report
↓
Serialize/report results
```

---

# Model Structure

Quadra models typically implement:

```cpp
template <typename Context>
void initialize(Context&)
{
}

template <typename T, typename Context>
T evaluate(const std::vector<T>& p, Context&) const
{
    return evaluate_impl<T>(p);
}
```

The model returns a scalar objective function value (negative log-likelihood).

---

# Workspace Lifecycle

The Laplace implicit workspace stores reusable inference quantities:

- random-effect Hessian
- mixed Hessian
- implicit derivatives
- sparse factorization cache
- optimized random effects

Example:

```cpp
const auto workspace =
    quadra::build_laplace_implicit_workspace(
        model,
        theta_hat,
        random_initial,
        parameters);
```

The workspace is reusable across:

- derived quantity uncertainty
- mode sensitivity analysis
- covariance propagation
- profiling workflows

---

# Fixed-Effect Covariance Estimation

Fixed-effect covariance estimation is performed from the Laplace objective:

```cpp
const auto covariance =
    quadra::estimate_fixed_effect_covariance(
        objective,
        theta_hat);
```

The resulting covariance matrix is used for:

- standard errors
- confidence intervals
- delta method propagation
- derived quantity uncertainty

---

# Profiled Derived Uncertainty

Quadra supports Laplace-aware derived uncertainty propagation using implicit
differentiation.

This accounts for:

```text
du / dtheta
```

which represents the sensitivity of optimized random effects to fixed effects.

Example:

```cpp
const auto report =
    quadra::compute_laplace_profiled_derived_report(
        quantities,
        theta_hat,
        workspace.u_hat_m,
        covariance.covariance_m,
        workspace);
```

The resulting report contains:

- estimates
- profiled Jacobian
- covariance matrix
- correlation matrix
- standard errors
- coefficients of variation

---

# CSV Export

Reports can be serialized directly to CSV:

```cpp
quadra::write_profiled_derived_report_csv(
    "report.csv",
    report);
```

Matrices can also be exported:

```cpp
quadra::write_matrix_csv(
    "covariance.csv",
    report.delta_m.covariance_m);
```

---

# Minimal Example Walkthrough

The canonical minimal V1 example is:

```text
examples/simple/random_intercept_model.cpp
```

This example demonstrates:

- parameter registration
- random effects
- Laplace approximation
- implicit differentiation
- covariance estimation
- profiled derived uncertainty
- CSV export
- timing diagnostics

---

# Catch-at-Age Example

The primary large-scale example is:

```text
examples/big/catch_at_age_laplace.cpp
```

This example demonstrates:

- high-dimensional random effects
- sparse Hessian factorization
- cached Laplace workflows
- profiled derived uncertainty
- realistic fisheries stock assessment structure

---

# Design Philosophy

Quadra emphasizes:

- composable inference primitives
- explicit workflows
- reusable workspaces
- sparse-aware numerical methods
- modern C++ design
- scalable scientific computing

The framework is intended to support both:

- methodological research
- production statistical modeling systems

---

# Current V1 Capabilities

Current V1 capabilities include:

- Laplace approximation
- sparse Hessian methods
- exact AD gradients
- implicit differentiation
- profiled delta method
- covariance estimation
- reusable factorization caches
- report serialization
- timing instrumentation
- mixed-effects workflows
