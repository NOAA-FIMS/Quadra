# Covariance estimation

Quadra now includes a first-pass fixed-effect covariance estimator in:

```text
core/inference/covariance.hpp
```

The current implementation computes:

```text
Cov(theta_hat) ≈ H_laplace(theta_hat)^-1
```

where `H_laplace` is approximated by central finite differences of the optimized Laplace gradient.

This is intentionally a validation-friendly first implementation. Future versions should replace the finite-difference Hessian of the Laplace gradient with exact second derivatives or selected sparse inverse machinery.

## Usage

```cpp
quadra::LaplaceOptions lopts;
lopts.use_hutchinson_trace = false;

auto fit = quadra::optimize_lbfgs(model, params, lopts);

quadra::CovarianceOptions copts;
copts.fd_step = 1e-4;

auto cov = quadra::estimate_fixed_covariance(
    model,
    params,
    fit,
    lopts,
    copts);
```

Standard errors are stored back onto parameters by default:

```cpp
params.params[i].std_error
params.params[i].covariance_index
```

## Current scope

Implemented:
- fixed-effect covariance,
- storage of standard errors on `Parameter`,
- covariance subset controlled by `Parameter::estimate_covariance`.

Not implemented yet:
- random-effect covariance,
- derived quantity covariance,
- delta method,
- selected inverse exact sparse trace.
