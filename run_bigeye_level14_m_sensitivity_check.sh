#!/usr/bin/env bash
set -euo pipefail

./inspect_bigeye_level14_m_sensitivity.sh

echo
echo "== O3 build Bigeye Level 14 M sensitivity =="
mkdir -p build/examples
c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  examples/NMFS/pifsc_bigeye_tuna/level14_m_sensitivity/quadra/bigeye_level14_m_sensitivity.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level14_m_sensitivity/quadra/bigeye_adgraph_global.cpp \
  -o build/examples/pifsc_bigeye_level14_m_sensitivity_check

OUTDIR="examples/NMFS/pifsc_bigeye_tuna/level14_m_sensitivity/outputs"
SUMMARY="$OUTDIR/bigeye_level14_m_sensitivity_summary.csv"
mkdir -p "$OUTDIR"
echo "m,objective,grad_norm,converged,r0,fbar,q_purse_seine,rec_sd,rec_lag1,rec_prior_nll,rec_max_abs,age10_multiplier,fitted_plus_bio_share,equilibrium_plus_bio_share" > "$SUMMARY"

for M in 0.18 0.25 0.35 0.45; do
  echo
  echo "== Run Bigeye Level 14 M sensitivity: M=$M =="
  BIGEYE_FIXED_M="$M" ./build/examples/pifsc_bigeye_level14_m_sensitivity_check

  cp "$OUTDIR/bigeye_level14_fit_summary.csv" "$OUTDIR/bigeye_level14_fit_summary_M_${M}.csv"
  cp "$OUTDIR/bigeye_level14_recruitment_diagnostics.csv" "$OUTDIR/bigeye_level14_recruitment_diagnostics_M_${M}.csv"
  cp "$OUTDIR/bigeye_level14_initial_numbers_diagnostics.csv" "$OUTDIR/bigeye_level14_initial_numbers_diagnostics_M_${M}.csv"

  python3 - "$M" "$SUMMARY" <<'PY'
import csv
import sys
from pathlib import Path

m = sys.argv[1]
summary_path = Path(sys.argv[2])
outdir = Path("examples/NMFS/pifsc_bigeye_tuna/level14_m_sensitivity/outputs")

fit = {}
with (outdir / "bigeye_level14_fit_summary.csv").open() as f:
    for r in csv.DictReader(f):
        fit[r["field"]] = r["value"]

rec = {}
with (outdir / "bigeye_level14_recruitment_diagnostics.csv").open() as f:
    for r in csv.DictReader(f):
        if r.get("section") == "summary":
            rec[r["metric"]] = r["value"]

init = {}
age = {}
with (outdir / "bigeye_level14_initial_numbers_diagnostics.csv").open() as f:
    for r in csv.DictReader(f):
        if r["section"] == "summary":
            init[r["metric"]] = r["value"]
        elif r["section"] == "age":
            age.setdefault(r["target"], {})[r["metric"]] = r["value"]

row = [
    m,
    fit.get("objective", ""),
    fit.get("grad_norm", ""),
    fit.get("converged", ""),
    fit.get("r0", ""),
    fit.get("fbar", ""),
    fit.get("q_purse_seine", ""),
    rec.get("sd", ""),
    rec.get("lag1_correlation", ""),
    rec.get("total_prior_nll", ""),
    rec.get("max_abs", ""),
    age.get("age_10", {}).get("multiplier", ""),
    init.get("fitted_plus_bio_share", ""),
    init.get("equilibrium_plus_bio_share", ""),
]
with summary_path.open("a", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(row)
PY
done

echo
echo "== Level 14 M sensitivity summary =="
cat "$SUMMARY"
