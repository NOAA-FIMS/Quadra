# Quadra performance framework

This directory defines the shared contract for Quadra correctness and
performance studies. Benchmarks should describe the model form they exercise,
emit one JSON object per case, and keep correctness gates separate from noisy
wall-clock regression gates.

## Model-form catalog

Every Laplace benchmark should cover or explicitly exclude these forms:

| Form | Required structural assertion | Typical backend |
|---|---|---|
| diagonal | bandwidth 0 | diagonal or banded |
| tridiagonal | bandwidth 1 | tridiagonal |
| banded | bounded bandwidth | banded |
| block diagonal | disconnected blocks | block/sparse |
| arrowhead | dense border with sparse body | sparse |
| general sparse | low fill ratio | sparse LDLT |
| nearly dense | high fill ratio | dense LDLT |
| dense | full pattern | dense LDLT |
| parameter-dependent | replay validation/rebuild | discovered backend |

`result.schema.json` is the stable artifact contract. Raw samples are retained;
summaries and regression decisions must be derived from them.

## Measurement rules

- Report cold recording and warm replay independently.
- Report objective, gradient, and Hessian agreement before performance.
- Record selected backend, dimensions, nonzeros, bandwidth, and tape rebuilds.
- Use fresh processes for peak RSS.
- Use medians of repeated samples for timing comparisons.
- Do not enforce wall-clock thresholds on unpinned shared CI runners.
- Ordinary CI may gate deterministic metrics such as backend choice, graph
  size, active vertices, allocations, and tape rebuild count.

Validate an artifact with:

```sh
python3 benchmarks/framework/validate_results.py results.jsonl
```

Build both catalogs with `make benchmark-model-forms`. The executables may
also be run independently so their JSONL output can be retained:

```sh
make benchmarks/framework/model_form_backend_benchmark
./benchmarks/framework/model_form_backend_benchmark 128 10 > backend.jsonl

make benchmarks/framework/laplace_model_catalog_benchmark
./benchmarks/framework/laplace_model_catalog_benchmark 128 10 > laplace.jsonl
```

The backend catalog isolates discovery and numerical factorization. The full
Laplace catalog additionally exercises public-API recording, Newton mode
solving, Hessian discovery, log determinants, warm replay, and tape stability.
Run every form in a fresh process when collecting attributable peak RSS:

```sh
python3 benchmarks/framework/run_laplace_catalog.py \
  --dimension 128 --repetitions 10 --output laplace-fresh.jsonl
```
