#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-20}"
LENGTHS="${2:-25,50,100,250,500,1000}"

STAMP="$(date +%Y%m%d_%H%M%S)"
OUTDIR="benchmarks/quadra_vs_tmb_runtime_suite_${STAMP}"
mkdir -p "$OUTDIR"

SUMMARY="${OUTDIR}/summary.tsv"
printf "source\tn\tobjective\tavg_ms\trelative_to_tmb\n" > "$SUMMARY"

run_capture() {
  local label="$1"
  local script="$2"
  local outfile="$3"

  echo
  echo "== ${label} =="

  if [[ -x "./${script}" ]]; then
    "./${script}" "$REPS" "$LENGTHS" | tee "${OUTDIR}/${outfile}"
  else
    echo "SKIP: ${script} not found or not executable" | tee "${OUTDIR}/${outfile}"
  fi
}

extract_tmb() {
  local file="$1"
  awk '
    /^[[:space:]]*[0-9]+[[:space:]]/ {
      print $1 "\t" $2 "\t" $3
    }
  ' "$file"
}

extract_laplace_eval_warm() {
  local file="$1"
  awk '
    /^[[:space:]]*[0-9]+[[:space:]]/ {
      # n objective lbfgs cold warm speedup cold_it warm_it grad diff
      print $1 "\t" $2 "\t" $5
    }
  ' "$file"
}

run_capture \
  "TMB scaled fixed-theta" \
  "run_tmb_scaled_fixed_theta_benchmark.sh" \
  "tmb_scaled_fixed_theta.txt"

run_capture \
  "Quadra LaplaceEvaluator warm-start runtime" \
  "run_state_space_laplace_evaluator_benchmark.sh" \
  "quadra_laplace_evaluator.txt"


TMB_TMP="${OUTDIR}/tmb_table.tsv"
Q_TMP="${OUTDIR}/quadra_laplace_evaluator_table.tsv"

: > "$TMB_TMP"
: > "$Q_TMP"

if [[ -f "${OUTDIR}/tmb_scaled_fixed_theta.txt" ]]; then
  extract_tmb "${OUTDIR}/tmb_scaled_fixed_theta.txt" > "$TMB_TMP" || true
fi

if [[ -f "${OUTDIR}/quadra_laplace_evaluator.txt" ]]; then
  extract_laplace_eval_warm "${OUTDIR}/quadra_laplace_evaluator.txt" > "$Q_TMP" || true
fi

# Write TMB rows first.
awk '
  NF == 3 {
    print "TMB_fixed_theta\t" $1 "\t" $2 "\t" $3 "\t1.0"
  }
' "$TMB_TMP" >> "$SUMMARY"

# Join Quadra rows to TMB by n.
awk '
  BEGIN {
    FS = OFS = "\t"
  }
  FNR == NR {
    tmb[$1] = $3
    next
  }
  NF == 3 {
    n=$1
    obj=$2
    ms=$3
    rel = (n in tmb && ms > 0.0) ? tmb[n] / ms : 0.0
    print "Quadra_LaplaceEvaluator_warm", n, obj, ms, rel
  }
' "$TMB_TMP" "$Q_TMP" >> "$SUMMARY"

cat > "${OUTDIR}/README.txt" <<README
Quadra vs TMB clean runtime suite
===========================

Generated: ${STAMP}
Reps:      ${REPS}
Lengths:   ${LENGTHS}

Files:
  tmb_scaled_fixed_theta.txt
  quadra_laplace_evaluator.txt
  summary.tsv

summary.tsv columns:
  source
  n
  objective
  avg_ms
  relative_to_tmb

relative_to_tmb:
  TMB rows are 1.0
  Quadra rows are TMB_ms / Quadra_ms
README

echo
echo "Quadra vs TMB suite complete."
echo "Output directory:"
echo "  ${OUTDIR}"
echo
echo "Summary:"
echo "  ${SUMMARY}"
echo
column -t -s $'\t' "$SUMMARY" || cat "$SUMMARY"

