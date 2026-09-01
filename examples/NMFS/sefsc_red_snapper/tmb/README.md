# TMB Reference Implementation

This directory contains the optional TMB comparison model for the synthetic
Red Snapper example. It mirrors the Quadra model for numerical validation; it
is not an official assessment or a standalone scientific reference.

From the repository root, run both implementations and produce the comparison
CSV with:

```bash
./examples/NMFS/sefsc_red_snapper/run_quadra_vs_tmb_comparison.sh
```

Requirements are `Rscript`, the R package `TMB`, Python 3, and a working R C++
toolchain. TMB compilation creates a platform-specific shared library beside
`red_snapper_tmb.cpp`; generated comparison products are written to `../outputs/`.
