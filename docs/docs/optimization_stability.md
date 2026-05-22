# Optimization Stability in Quadra

## Overview

Quadra uses a nested optimization structure for Laplace approximation:

- Outer optimization over fixed effects using L-BFGS
- Inner optimization over random effects using sparse Newton solves
- Exact or approximate Laplace objective evaluation at the random-effect mode

Large mixed-effects models can encounter unstable trial evaluations during optimization. The goal of the recent stabilization work was not to eliminate all failed trial evaluations, but rather to:

1. Recover safely from them
2. Preserve optimizer progress
3. Avoid false convergence
4. Maintain reproducibility across platforms (macOS/Linux)

---

# Adaptive Newton Damping

## Problem

Sparse Newton solves may fail when the Hessian is:

- indefinite
- poorly conditioned
- numerically singular
- outside the local quadratic region

This commonly appeared as:

Sparse Hessian factorization failed

during random-effect mode solving.

## Solution

Quadra now applies adaptive diagonal damping during sparse Newton solves.

(H + λI)s = -g

with progressively increasing damping values until:

- factorization succeeds
- the step is finite
- the solve becomes numerically stable

## Benefits

- improved robustness
- reduced catastrophic failures
- preserved convergence behavior
- improved Linux/macOS consistency

---

# Recoverable Laplace Trial Failures

Not every failed trial evaluation represents optimizer failure.

Quadra now catches recoverable failures during:

- random-effect mode solves
- Laplace objective evaluation

and returns a large penalty objective value instead of aborting optimization.

This allows L-BFGS to reject unstable trial steps and continue searching.

---

# Directional Penalty Gradients

Earlier penalty paths returned:

- large objective
- zero gradient

This can falsely signal convergence.

Penalty paths now return directional gradients proportional to the current parameter vector, encouraging movement away from unstable regions.

---

# Why Penalty Counts Are Allowed

The convergence contract intentionally allows some penalty-path evaluations.

Example:

NOTE: Laplace penalty path occurred during trial evaluations: 5

This is acceptable because the final solution:

- converged successfully
- produced finite gradients
- passed objective consistency checks

---

# Reverted Step-Norm Safeguard

An experimental Newton step clipping safeguard increased instability on Linux CI and was reverted.

Lesson:

Numerically conservative safeguards can still destabilize quasi-Newton trajectories.

---

# Future Work

Potential future directions:

- trust-region Newton methods
- ratio-based acceptance tests
- Hessian diagnostics
- adaptive line-search heuristics

---

# Key Principle

Quadra prioritizes robust recovery and truthful convergence reporting over forcing every intermediate evaluation to succeed.
