# Correlated random-effect benchmark

This benchmark tests sparse structure beyond independent random effects.

Observation model:

```text
y_ig ~ N(mu + u_g, 1)
```

AR1-style random-effect prior:

```text
0.5 * lambda0 * u_0^2
+ 0.5 * lambda_diff * sum_{g=1}^{G-1} (u_g - rho*u_{g-1})^2
```

Expected random-effect Hessian pattern:

```text
tridiagonal
nnz = 3G - 2
```

Run:

```bash
make tests/test_correlated_random_intercept_hessian
./tests/test_correlated_random_intercept_hessian

make benchmark-correlated-random
```
