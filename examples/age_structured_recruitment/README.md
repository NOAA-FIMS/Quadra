# Age-structured recruitment deviation benchmark

This example tests whether the Quadra structure-aware advantage carries beyond
surplus production.

The model has:

```text
numbers at age
constant natural mortality
logistic survey selectivity
survey index likelihood
annual recruitment deviations
```

First version:

```text
x[y] ~ Normal(0, sigma_R)
```

so the recruitment-deviation Hessian is diagonal. The next version should use
RW1/AR1 recruitment deviations, yielding a tridiagonal Hessian.

Run:

```bash
bash examples/age_structured_recruitment/run_quadra_vs_tmb_age_structured_recruitment_benchmark.sh \
  10 25,50,100,250,500,1000 10
```
