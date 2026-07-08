#!/usr/bin/env bash
set -euo pipefail

SRC="examples/NMFS/pifsc_bigeye_tuna/level22_initial_number_regularization_scan"
CPP="$SRC/quadra/bigeye_level22_initial_number_regularization_scan.cpp"
BIN="/tmp/bigeye_level22_initial_number_regularization_scan"

if [[ ! -f "$CPP" ]]; then
  echo "ERROR: missing $CPP"
  exit 1
fi

OUT="$SRC/outputs"
mkdir -p "$OUT"

SUMMARY="$OUT/bigeye_level22_initial_number_sigma_scan_summary.csv"
echo "sigma_init,objective,grad_norm,converged,longline_mean_abs,longline_max_abs,purse_seine_mean_abs,purse_seine_max_abs,initial_numbers_prior_nll,longline_selectivity_prior_nll,purse_seine_selectivity_prior_nll,recruitment_prior_nll" > "$SUMMARY"

compile() {
  echo
  echo "== O3 build Bigeye Level 22 initial-number regularization scan =="
  clang++ -std=c++17 -O3 -I. -Icore -Icore/eigen -Iexternal/eigen \
    "$CPP" "$SRC/quadra/bigeye_adgraph_global.cpp" \
    -o "$BIN"
}

extract_field() {
  local file="$1"
  local key="$2"
  awk -F, -v k="$key" '$1==k {print $2; exit}' "$file"
}

extract_residual_metric() {
  local file="$1"
  local fleet="$2"
  local col="$3"
  awk -F, -v fleet="$fleet" -v col="$col" '
    $1==fleet {
      if (col=="mean") print $2;
      else if (col=="max") print $3;
      exit
    }' "$file"
}

compile

for sig in 0.35 0.50 0.75 1.00 1.25; do
  echo
  echo "== Run Bigeye Level 22 initial-number sigma=${sig} =="
  LEVEL22_SIGMA_INIT_DEV="$sig" "$BIN"

  FIT="$OUT/bigeye_level22_fit_summary.csv"
  RES="$OUT/bigeye_level22_age_comp_residual_diagnostics.csv"
  SAN="$OUT/bigeye_level22_parameter_sanity_diagnostics.csv"

  objective="$(extract_field "$FIT" objective)"
  grad_norm="$(extract_field "$FIT" grad_norm)"
  converged="$(extract_field "$FIT" converged)"

  ll_mean="$(extract_residual_metric "$RES" longline mean)"
  ll_max="$(extract_residual_metric "$RES" longline max)"
  ps_mean="$(extract_residual_metric "$RES" purse_seine mean)"
  ps_max="$(extract_residual_metric "$RES" purse_seine max)"

  init_prior="$(extract_field "$SAN" initial_numbers_prior_nll)"
  ll_prior="$(extract_field "$SAN" longline_selectivity_prior_nll)"
  ps_prior="$(extract_field "$SAN" purse_seine_selectivity_prior_nll)"
  rec_prior="$(extract_field "$SAN" recruitment_prior_nll)"

  echo "$sig,$objective,$grad_norm,$converged,$ll_mean,$ll_max,$ps_mean,$ps_max,$init_prior,$ll_prior,$ps_prior,$rec_prior" >> "$SUMMARY"

  cp "$FIT" "$OUT/bigeye_level22_fit_summary_sigma_${sig}.csv"
  cp "$RES" "$OUT/bigeye_level22_age_comp_residual_diagnostics_sigma_${sig}.csv"
  cp "$SAN" "$OUT/bigeye_level22_parameter_sanity_diagnostics_sigma_${sig}.csv"
done

echo
echo "== Level 22 initial-number regularization scan summary =="
cat "$SUMMARY"

echo
echo "== Best rows sorted by objective =="
{ head -1 "$SUMMARY"; tail -n +2 "$SUMMARY" | sort -t, -k2,2n; } | head -8
