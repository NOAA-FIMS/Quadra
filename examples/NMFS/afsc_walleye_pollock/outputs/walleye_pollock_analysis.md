# Synthetic Walleye Pollock Functional Analysis

Synthetic and public-data-safe. Not an official assessment.

## Executive Summary

- **Overall status:** `HEALTHY`.
- **Confidence:** `HIGH`.
- **Optimization quality:** `EXCELLENT`.
- **Uncertainty structure:** `LOCAL`.
- **Optimization:** converged = `yes`, gradient norm = `5.19827398754203e-06`.
- **Curvature health:** positive definite = `yes`, condition number = `10.9840248198155`.
- **Latent structure:** `20` random effects were estimated.
- **Symbolic vs numerical structure:** structural density = `0.91`, but 95% of curvature is retained by `58` entries.
- **Spectral complexity:** entropy effective rank = `16.3107626636197`, with 90% curvature requiring `14` eigen-directions.

## Model Health Assessment

| Check | Status | Evidence |
|---|---:|---|
| Optimization | `PASS` | converged = `yes` |
| Gradient quality | `PASS` | gradient norm = `5.19827398754203e-06` |
| Curvature | `PASS` | positive definite = `yes` |
| Conditioning | `EXCELLENT` | condition number = `10.9840248198155` |
| Overall status | `HEALTHY` | rule-based v1 diagnostic |
| Confidence | `HIGH` | based on convergence, gradient, PD status, and conditioning |

**Interpretation:** the rule-based health check is intentionally simple. It flags obvious numerical issues quickly, but it does not replace scientific review or model-specific diagnostics.

## Model Complexity

| Quantity | Value |
|---|---:|
| Fixed effects | `2` |
| Random effects | `20` |
| Total estimated quantities | `22` |
| Structural nonzeros | `364` |
| Structural density | `0.91` |
| Entries for 95% curvature | `58` |
| Effective bandwidth for 95% curvature | `1` |
| 95% curvature compression | `6.27586x` |

## Optimization

- Quality: `EXCELLENT`
- Objective value: `-14.3868675854755`
- Gradient norm: `5.19827398754203e-06`
- Converged: `yes`
- Max gradient parameter: `log_r0`

## Curvature

- Positive definite: `yes`
- Condition number: `10.9840248198155`
- Minimum eigenvalue: `10.2183884027573`
- Maximum eigenvalue: `112.239031834401`

## Spectral Structure

- Largest eigenvalue share: `0.095050325395682`
- Entropy effective rank: `16.3107626636197`
- Eigenvectors needed for 90% curvature: `14`
- Eigenvectors needed for 95% curvature: `16`

**Interpretation:** curvature is distributed across many latent-state directions rather than being dominated by one or two modes. That is a good sign for numerical stability.

## Effective Structure

- Structural density: `0.91`
- Structural nonzeros: `364`
- Entries for 95% curvature: `58`
- Effective bandwidth for 95% curvature: `1`
- 95% curvature compression: `6.27586x`

**Interpretation:** symbolic density alone overstates practical complexity. The detailed Laplace report below shows that large amounts of curvature can be retained with far fewer entries or a narrow effective bandwidth.

## Correlation Graph

- Classification: `LOCAL`
- Average degree: `1.5`
- Maximum degree: `2`
- Connected components: `5`
- Largest component size: `14`
- Graph diameter: `13`

**Interpretation:** a LOCAL graph means the strongest uncertainty relationships are neighborhood-like rather than globally tangled.

## Latent State Summary

- Count: `20`
- Mean: `0.0824010336165072`
- Standard deviation: `0.0500130055529492`

## Key Takeaway

This report demonstrates why Quadra's functional analysis diagnostics are useful: a model can look dense from a symbolic Hessian pattern, while numerical curvature, graph structure, and effective bandwidth reveal a simpler local-dependence structure.

## Full Laplace Structure Report

```text
Laplace Structure Report
========================

Random effects:              20
Matrix size:                 20 x 20
Total entries:               400
Structural nonzeros:         364 / 400 (91%)
Nonzero tolerance:           1e-08
Max |H_ij|:                  61.5960537686533
Positive definite:           yes
Min eigenvalue:              10.2183884027573
Max eigenvalue:              112.239031834401
Condition number:            10.9840248198155

Effective sparsity
------------------
curvature_retained,entries_required,entry_share,compression_vs_structural
90%,54,0.135,6.74074074074074
95%,58,0.145,6.27586206896552
97%,100,0.25,3.64
98%,133,0.3325,2.73684210526316
99%,183,0.4575,1.98907103825137
99.5%,224,0.56,1.625
99.9%,293,0.7325,1.24232081911263
100%,364,0.91,1

Effective bandwidth
-------------------
curvature_retained,bandwidth,entry_count_if_banded,entry_share_if_banded
90%,1,58,0.145
95%,1,58,0.145
97%,2,94,0.235
98%,3,128,0.32
99%,5,190,0.475
99.5%,7,244,0.61
99.9%,10,310,0.775
100%,19,400,1

Interpretation
--------------
This report measures numerical curvature concentration, not only symbolic sparsity.
A dense structural Hessian can still be effectively sparse if most curvature is carried by relatively few entries or bands.
```
