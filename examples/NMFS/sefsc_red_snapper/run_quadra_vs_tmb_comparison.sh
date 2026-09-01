#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$repo_root"

command -v Rscript >/dev/null 2>&1 || {
  echo "error: Rscript is required for the TMB comparison" >&2
  exit 1
}
command -v python3 >/dev/null 2>&1 || {
  echo "error: python3 is required for the comparison report" >&2
  exit 1
}

./examples/NMFS/sefsc_red_snapper/run_red_snapper_quadra_fit.sh
Rscript examples/NMFS/sefsc_red_snapper/tmb/run_red_snapper_tmb_fit.R
python3 examples/NMFS/sefsc_red_snapper/compare_quadra_tmb_fit.py
cat examples/NMFS/sefsc_red_snapper/outputs/quadra_vs_tmb_fit_comparison.csv
