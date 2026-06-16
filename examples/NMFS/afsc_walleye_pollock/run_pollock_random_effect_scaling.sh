#!/usr/bin/env bash
set -euo pipefail

echo "Assessment-scale synthetic diagnostics: using QUADRA_LBFGS_GRAD_TOL=1.0e-2"
echo "This scaling run is intended to exercise assessment-like random-effect behavior, not strict optimizer validation."
echo

mkdir -p build/examples examples/NMFS/afsc_walleye_pollock/outputs

CPP="examples/NMFS/afsc_walleye_pollock/quadra/walleye_pollock.cpp"
GLOB="examples/NMFS/afsc_walleye_pollock/quadra/walleye_pollock_adgraph_global.cpp"
OUT="examples/NMFS/afsc_walleye_pollock/outputs/random_effect_scaling_summary.csv"

echo "random_effects,exit_code,objective,grad_norm,converged,max_grad_param,max_grad_value,max_abs_grad,message" > "$OUT"

run_case() {
  local n="$1"
  local exe="build/examples/afsc_walleye_pollock_re_${n}"
  local log="examples/NMFS/afsc_walleye_pollock/outputs/random_effect_scaling_${n}.log"

  echo
  echo "== Random recruitment count: ${n} =="

  rm -f examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_fit_summary.csv
  rm -f examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_recruitment_deviations.csv
  rm -f examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_fixed_gradient_diagnostics.csv
  rm -f examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_fixed_parameter_estimates.csv
  rm -f examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_fixed_hessian_diagnostics.csv
  rm -f examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_fixed_hessian_matrix.csv
  rm -f examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_fixed_hessian_diagnostics.csv
  rm -f examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_fixed_hessian_matrix.csv
  rm -f examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_huu_diagnostics.csv

  if [[ "$n" == "0" ]]; then
    clang++ -std=c++17 -g -I"external/eigen/"       -DWALLEYE_POLLOCK_HUU_DIAGNOSTICS       -DQUADRA_LBFGS_GRAD_TOL=1.0e-2       "$CPP" "$GLOB" -o "$exe"
  else
    clang++ -std=c++17 -g -I"external/eigen/"       -DWALLEYE_POLLOCK_HUU_DIAGNOSTICS       -DQUADRA_LBFGS_GRAD_TOL=1.0e-2       -DWALLEYE_POLLOCK_RANDOM_RECRUITMENT_COUNT="$n"       "$CPP" "$GLOB" -o "$exe"
  fi

  set +e
  "$exe" > "$log" 2>&1
  local code="$?"
  set -e

  cat "$log"

  local huu_diag="examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_huu_diagnostics.csv"
  if [[ -f "$huu_diag" ]]; then
    echo
    echo "Huu diagnostics:"
    cat "$huu_diag"
  fi

  local param_diag="examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_fixed_parameter_estimates.csv"
  if [[ -f "$param_diag" ]]; then
    echo
    echo "Fixed-parameter estimates:"
    cat "$param_diag"
  fi

  local grad_diag="examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_fixed_gradient_diagnostics.csv"
  if [[ -f "$grad_diag" ]]; then
    echo
    echo "Fixed-gradient diagnostics:"
    cat "$grad_diag"
  fi

  local fixed_hess_diag="examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_fixed_hessian_diagnostics.csv"
  if [[ -f "$fixed_hess_diag" ]]; then
    echo
    echo "Fixed-effect Hessian diagnostics:"
    cat "$fixed_hess_diag"
  fi

  local fixed_hess_diag="examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_fixed_hessian_diagnostics.csv"
  if [[ -f "$fixed_hess_diag" ]]; then
    echo
    echo "Fixed-effect Hessian diagnostics:"
    cat "$fixed_hess_diag"
  fi

  local summary="examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_fit_summary.csv"
  local objective="NA"
  local grad_norm="NA"
  local converged="no"
  local message="run_failed"
  local max_grad_param="NA"
  local max_grad_value="NA"
  local max_abs_grad="NA"

  if [[ -f "$grad_diag" ]]; then
    max_grad_param="$(awk -F, 'NR>1 {if ($3+0 > max) {max=$3+0; p=$1; g=$2}} END {if (p!="") print p; else print "NA"}' "$grad_diag")"
    max_grad_value="$(awk -F, 'NR>1 {if ($3+0 > max) {max=$3+0; p=$1; g=$2}} END {if (p!="") print g; else print "NA"}' "$grad_diag")"
    max_abs_grad="$(awk -F, 'NR>1 {if ($3+0 > max) max=$3+0} END {if (max!="") print max; else print "NA"}' "$grad_diag")"
  fi

  if [[ -f "$summary" ]]; then
    objective="$(awk -F, '$1=="objective"{print $2}' "$summary" | tail -1)"
    grad_norm="$(awk -F, '$1=="grad_norm"{print $2}' "$summary" | tail -1)"
    converged="$(awk -F, '$1=="converged"{print $2}' "$summary" | tail -1)"
    message="$(awk -F, '$1=="message"{print $2}' "$summary" | tail -1)"
  fi

  echo "${n},${code},${objective},${grad_norm},${converged},${max_grad_param},${max_grad_value},${max_abs_grad},${message}" >> "$OUT"
}

for n in 0 1 2 5 10 20; do
  run_case "$n"
done

echo
echo "Wrote scaling summary:"
echo "  $OUT"
cat "$OUT"
