# LaplaceEvaluator

This patch introduces the first conservative `LaplaceEvaluator`.

It caches:
- `ParameterPartition`
- warm-started `u_hat`
- `CachedLaplaceObjectiveState`
- sparse factorization symbolic analysis

It does not yet cache reusable AD tapes.

Run:

```bash
make tests/test_laplace_evaluator
./tests/test_laplace_evaluator

make benchmark-laplace-evaluator
```
