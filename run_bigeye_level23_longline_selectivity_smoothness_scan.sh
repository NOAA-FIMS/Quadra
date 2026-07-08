#!/usr/bin/env bash
set -euo pipefail

SRC="examples/NMFS/pifsc_bigeye_tuna/level23_longline_selectivity_smoothness_scan"
CPP="$SRC/quadra/bigeye_level23_longline_selectivity_smoothness_scan.cpp"
BIN="/tmp/bigeye_level23_longline_selectivity_smoothness_scan"

if [[ ! -f "$CPP" ]]; then
  echo "ERROR: missing $CPP"
  exit 1
fi

OUT="$SRC/outputs"
mkdir -p "$OUT"

SUMMARY="$OUT/bigeye_level23_longline_smoothness_scan_summary.csv"
echo "ll_sigma,smooth_lambda,objective,grad_norm,converged,longline_mean_abs,longline_max_abs,purse_seine_mean_abs,purse_seine_max_abs,longline_selectivity_prior_nll,initial_numbers_prior_nll,recruitment_prior_nll" > "$SUMMARY"

echo
echo "== O3 build Bigeye Level 23 longline selectivity smoothness scan =="
clang++ -std=c++17 -O3 -I. -Icore -Icore/eigen -Iexternal/eigen \
  "$CPP" "$SRC/quadra/bigeye_adgraph_global.cpp" \
  -o "$BIN"

extract_field() {
  local file="$1"
  local key="$2"
  awk -F, -v k="$key" '$1==k {print $2; exit}' "$file"
}

extract_residual_metric() {
  local file="$1"
  local fleet="$2"
  local col="$3"

  # Prefer the text report's Fleet summary block because it is stable:
  # fleet,mean_abs_residual,max_abs_residual,n
  awk -F, -v fleet="$fleet" -v col="$col" '
    /^Fleet summary/ { in_block=1; next }
    in_block && /^$/ { in_block=0 }
    in_block && $1==fleet {
      if (col=="mean") print $2;
      else if (col=="max") print $3;
      exit
    }
  ' "${file%.csv}.txt"
}

extract_named_csv_value() {
  local file="$1"
  local key="$2"

  # Parameter sanity writes a text block with lines like:
  # longline_selectivity_prior_nll,11.18
  awk -F, -v k="$key" '$1==k { print $2; exit }' "${file%.csv}.txt"
}

for smooth in 0 0.01 0.05 0.10 0.25 0.50 1.00; do
  echo
  echo "== Run Bigeye Level 23 ll_sigma=1.75 smooth_lambda=${smooth} =="
  LEVEL23_LL_SEL_SIGMA="1.75" LEVEL23_LL_SMOOTH_LAMBDA="$smooth" "$BIN"

  FIT="$OUT/bigeye_level23_fit_summary.csv"
  RES="$OUT/bigeye_level23_age_comp_residual_diagnostics.csv"
  SAN="$OUT/bigeye_level23_parameter_sanity_diagnostics.csv"

  objective="$(extract_field "$FIT" objective)"
  grad_norm="$(extract_field "$FIT" grad_norm)"
  converged="$(extract_field "$FIT" converged)"

  ll_mean="$(extract_residual_metric "$RES" longline mean)"
  ll_max="$(extract_residual_metric "$RES" longline max)"
  ps_mean="$(extract_residual_metric "$RES" purse_seine mean)"
  ps_max="$(extract_residual_metric "$RES" purse_seine max)"

  ll_prior="$(extract_named_csv_value "$SAN" longline_selectivity_prior_nll)"
  init_prior="$(extract_named_csv_value "$SAN" initial_numbers_prior_nll)"
  rec_prior="$(extract_named_csv_value "$SAN" recruitment_prior_nll)"

  echo "1.75,$smooth,$objective,$grad_norm,$converged,$ll_mean,$ll_max,$ps_mean,$ps_max,$ll_prior,$init_prior,$rec_prior" >> "$SUMMARY"

  cp "$FIT" "$OUT/bigeye_level23_fit_summary_smooth_${smooth}.csv"
  cp "$RES" "$OUT/bigeye_level23_age_comp_residual_diagnostics_smooth_${smooth}.csv"
  cp "$SAN" "$OUT/bigeye_level23_parameter_sanity_diagnostics_smooth_${smooth}.csv"

  if [[ -f "$OUT/bigeye_level23_longline_prediction_decomposition.csv" ]]; then
    cp "$OUT/bigeye_level23_longline_prediction_decomposition.csv" \
      "$OUT/bigeye_level23_longline_prediction_decomposition_smooth_${smooth}.csv"
  fi
  if [[ -f "$OUT/bigeye_level23_purse_seine_prediction_decomposition.csv" ]]; then
    cp "$OUT/bigeye_level23_purse_seine_prediction_decomposition.csv" \
      "$OUT/bigeye_level23_purse_seine_prediction_decomposition_smooth_${smooth}.csv"
  fi
done

echo
echo "== Level 23 longline smoothness scan summary =="
cat "$SUMMARY"

echo
echo "== Best rows sorted by objective =="
{ head -1 "$SUMMARY"; tail -n +2 "$SUMMARY" | sort -t, -k3,3n; } | head -10
