#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-20}"
LENGTHS="${2:-1000,2000,5000,10000,20000}"

STAMP="$(date +%Y%m%d_%H%M%S)"
OUTDIR="benchmarks/state_space_laplace_memory_${STAMP}"
mkdir -p "$OUTDIR"

BIN="build/examples/benchmark_latent_tridiagonal_laplace_evaluator"
SRC="examples/state_space_surplus_production/benchmark_latent_tridiagonal_laplace_evaluator.cpp"

if [[ ! -f "$SRC" ]]; then
  echo "ERROR: missing $SRC"
  exit 1
fi

mkdir -p build/examples

echo "== Building LaplaceEvaluator benchmark =="
set -x
c++ -std=c++17 -O2 -DNDEBUG -g \
  -Iexternal/Eigen \
  -Iexternal/LBFGSpp/include \
  -I. \
  -Iexamples/state_space_surplus_production \
  -Iexamples/surplus_production \
  "$SRC" \
  -o "$BIN"
set +x

echo
echo "== Running memory benchmark =="
echo "reps = ${REPS}"
echo "lengths = ${LENGTHS}"
echo

TIME_RAW="${OUTDIR}/raw_time.txt"
BENCH_OUT="${OUTDIR}/benchmark_output.txt"
SUMMARY="${OUTDIR}/summary.txt"

if /usr/bin/time -l true >/dev/null 2>&1; then
  TIME_MODE="macos_time_l"
  { /usr/bin/time -l "$BIN" "$REPS" "$LENGTHS" > "$BENCH_OUT"; } 2> "$TIME_RAW"
elif /usr/bin/time -v true >/dev/null 2>&1; then
  TIME_MODE="linux_time_v"
  { /usr/bin/time -v "$BIN" "$REPS" "$LENGTHS" > "$BENCH_OUT"; } 2> "$TIME_RAW"
else
  TIME_MODE="none"
  "$BIN" "$REPS" "$LENGTHS" | tee "$BENCH_OUT"
  echo "No supported /usr/bin/time RSS mode found." > "$TIME_RAW"
fi

cat "$BENCH_OUT"

MAX_RSS_BYTES=""
MAX_RSS_KB=""

if [[ "$TIME_MODE" == "macos_time_l" ]]; then
  MAX_RSS_BYTES="$(awk '/maximum resident set size/ {print $1; exit}' "$TIME_RAW" || true)"
  if [[ -n "$MAX_RSS_BYTES" ]]; then
    MAX_RSS_KB="$(awk -v b="$MAX_RSS_BYTES" 'BEGIN { printf "%.3f", b / 1024.0 }')"
  fi
elif [[ "$TIME_MODE" == "linux_time_v" ]]; then
  MAX_RSS_KB="$(awk -F: '/Maximum resident set size/ {gsub(/[ \t]/, "", $2); print $2; exit}' "$TIME_RAW" || true)"
  if [[ -n "$MAX_RSS_KB" ]]; then
    MAX_RSS_BYTES="$(awk -v kb="$MAX_RSS_KB" 'BEGIN { printf "%.0f", kb * 1024.0 }')"
  fi
fi

{
  echo "State-space LaplaceEvaluator memory benchmark"
  echo "============================================"
  echo
  echo "timestamp: ${STAMP}"
  echo "reps: ${REPS}"
  echo "lengths: ${LENGTHS}"
  echo "time_mode: ${TIME_MODE}"
  echo
  if [[ -n "${MAX_RSS_KB}" ]]; then
    echo "max_rss_kb: ${MAX_RSS_KB}"
    awk -v kb="$MAX_RSS_KB" 'BEGIN { printf "max_rss_mb: %.3f\n", kb / 1024.0 }'
  else
    echo "max_rss_kb: unavailable"
    echo "max_rss_mb: unavailable"
  fi
  echo
  echo "benchmark_output: ${BENCH_OUT}"
  echo "time_raw: ${TIME_RAW}"
} | tee "$SUMMARY"

echo
echo "Memory benchmark complete."
echo "Output directory:"
echo "  ${OUTDIR}"
