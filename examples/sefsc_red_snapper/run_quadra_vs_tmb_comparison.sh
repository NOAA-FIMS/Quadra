#!/usr/bin/env bash
set -euo pipefail

./examples/sefsc_red_snapper/run_red_snapper_quadra_fit.sh
Rscript examples/sefsc_red_snapper/tmb/run_red_snapper_tmb_fit.R
python3 examples/sefsc_red_snapper/compare_quadra_tmb_fit.py
cat examples/sefsc_red_snapper/outputs/quadra_vs_tmb_fit_comparison.csv
