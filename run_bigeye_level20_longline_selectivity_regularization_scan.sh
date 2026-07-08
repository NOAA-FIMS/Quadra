#!/usr/bin/env bash
set -euo pipefail

SRC="examples/NMFS/pifsc_bigeye_tuna/level20_longline_selectivity_regularization_scan/quadra/bigeye_level20_longline_selectivity_regularization_scan.cpp"
ADG="examples/NMFS/pifsc_bigeye_tuna/level20_longline_selectivity_regularization_scan/quadra/bigeye_adgraph_global.cpp"
OUTDIR="examples/NMFS/pifsc_bigeye_tuna/level20_longline_selectivity_regularization_scan/outputs"
SUMMARY="$OUTDIR/bigeye_level20_ll_selectivity_sigma_scan_summary.csv"
TXT="$OUTDIR/bigeye_level20_ll_selectivity_sigma_scan_summary.txt"

mkdir -p build/examples "$OUTDIR"

echo "sigma,objective,grad_norm,converged,longline_mean_abs,longline_max_abs,purse_seine_mean_abs,purse_seine_max_abs,longline_selectivity_prior_nll,initial_numbers_prior_nll,recruitment_prior_nll" > "$SUMMARY"
{
  echo "Level 20 Longline Selectivity Regularization Scan"
  echo "================================================"
  echo
  echo "sigma,objective,grad_norm,converged,longline_mean_abs,longline_max_abs,purse_seine_mean_abs,purse_seine_max_abs,longline_selectivity_prior_nll,initial_numbers_prior_nll,recruitment_prior_nll"
} > "$TXT"

for sigma in 0.50 0.75 1.00 1.25; do
  tag="${sigma/./p}"
  bin="build/examples/pifsc_bigeye_level20_ll_sigma_${tag}"

  echo
  echo "== O3 build Bigeye Level 20 longline selectivity sigma=${sigma} =="
  c++ -std=c++17 -O3 -flto \
    -DBIGEYE_LL_SEL_SIGMA=${sigma} \
    -I. -Iexternal/eigen \
    "$SRC" "$ADG" \
    -o "$bin"

  echo
  echo "== Run Bigeye Level 20 longline selectivity sigma=${sigma} =="
  "$bin" | tee "$OUTDIR/bigeye_level20_sigma_${tag}_run.log"

  cp "$OUTDIR/bigeye_level20_fit_summary.csv" "$OUTDIR/bigeye_level20_sigma_${tag}_fit_summary.csv"
  cp "$OUTDIR/bigeye_level20_age_comp_residual_diagnostics.txt" "$OUTDIR/bigeye_level20_sigma_${tag}_age_comp_residual_diagnostics.txt"
  cp "$OUTDIR/bigeye_level20_parameter_sanity_diagnostics.txt" "$OUTDIR/bigeye_level20_sigma_${tag}_parameter_sanity_diagnostics.txt"
  cp "$OUTDIR/bigeye_level20_parameter_sanity_diagnostics.csv" "$OUTDIR/bigeye_level20_sigma_${tag}_parameter_sanity_diagnostics.csv"

  python3 - "$sigma" "$SUMMARY" "$TXT" <<'PY'
from pathlib import Path
import csv
import sys

sigma, summary_path, txt_path = sys.argv[1], Path(sys.argv[2]), Path(sys.argv[3])
outdir = Path("examples/NMFS/pifsc_bigeye_tuna/level20_longline_selectivity_regularization_scan/outputs")

def read_summary(path):
    d = {}
    with open(path) as f:
        for row in csv.reader(f):
            if len(row) >= 2:
                d[row[0]] = row[1]
    return d

def fleet_summary(path):
    d = {}
    lines = Path(path).read_text().splitlines()
    in_block = False
    for line in lines:
        if line.startswith("fleet,mean_abs_residual"):
            in_block = True
            continue
        if in_block:
            if not line.strip():
                break
            parts = line.split(",")
            if len(parts) >= 4:
                d[parts[0]] = {"mean_abs": parts[1], "max_abs": parts[2], "n": parts[3]}
    return d

def prior_blocks(path):
    d = {}
    lines = Path(path).read_text().splitlines()
    in_block = False
    for line in lines:
        if line.startswith("Prior penalty by block"):
            in_block = True
            continue
        if in_block:
            if not line.strip() or line.startswith("Initial-number"):
                break
            if "," in line:
                k, v = line.split(",", 1)
                d[k] = v
    return d

fit = read_summary(outdir / "bigeye_level20_fit_summary.csv")
fleet = fleet_summary(outdir / "bigeye_level20_age_comp_residual_diagnostics.txt")
prior = prior_blocks(outdir / "bigeye_level20_parameter_sanity_diagnostics.txt")

row = [
    sigma,
    fit.get("objective", ""),
    fit.get("grad_norm", ""),
    fit.get("converged", ""),
    fleet.get("longline", {}).get("mean_abs", ""),
    fleet.get("longline", {}).get("max_abs", ""),
    fleet.get("purse_seine", {}).get("mean_abs", ""),
    fleet.get("purse_seine", {}).get("max_abs", ""),
    prior.get("longline_selectivity_prior_nll", ""),
    prior.get("initial_numbers_prior_nll", ""),
    prior.get("recruitment_prior_nll", ""),
]

with summary_path.open("a") as f:
    f.write(",".join(row) + "\n")
with txt_path.open("a") as f:
    f.write(",".join(row) + "\n")
PY

done

echo
echo "== Level 20 longline selectivity regularization scan summary =="
cat "$SUMMARY"

echo
echo "== Best rows sorted by objective =="
python3 - <<'PY'
from pathlib import Path
import csv

path = Path("examples/NMFS/pifsc_bigeye_tuna/level20_longline_selectivity_regularization_scan/outputs/bigeye_level20_ll_selectivity_sigma_scan_summary.csv")
rows = list(csv.DictReader(path.open()))
rows.sort(key=lambda r: float(r["objective"]))
for r in rows:
    print(
        f"sigma={r['sigma']}, objective={r['objective']}, "
        f"ll_mean_abs={r['longline_mean_abs']}, ll_max_abs={r['longline_max_abs']}, "
        f"ll_prior={r['longline_selectivity_prior_nll']}, "
        f"init_prior={r['initial_numbers_prior_nll']}, rec_prior={r['recruitment_prior_nll']}"
    )
PY
