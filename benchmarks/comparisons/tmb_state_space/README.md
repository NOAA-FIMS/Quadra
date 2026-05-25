# TMB State-Space Comparison

This benchmark compares Quadra and TMB on a Gaussian random-walk state-space model.

The comparison is intended to stress:

- latent state scaling
- sparse Hessian structure
- Laplace approximation cost
- derivative evaluation cost
- factorization behavior

The TMB benchmark is optional and requires R with TMB installed.
