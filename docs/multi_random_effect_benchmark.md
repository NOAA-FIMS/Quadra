# Multi-random-effect benchmark

This benchmark moves beyond the single-random-effect case.

Model:

```text
y_ij ~ N(mu + u_j, 1)
u_j  ~ N(0, 1)
```

It scales the number of random effects/groups `G` while holding observations per
group fixed.

This tests:
- fixed/random partitioning with many random effects
- sparse Hessian wrt random effects
- random-effect Newton solves
- Laplace objective/logdet
- exact-envelope LBFGS optimization

Run:

```bash
make tests/test_multi_random_intercept_laplace
./tests/test_multi_random_intercept_laplace

make benchmark-multi-random
```
