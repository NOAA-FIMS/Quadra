#!/usr/bin/env bash
set -euo pipefail

reps="${1:-10}"
lengths="${2:-25,50,100,250,500,1000}"
artifact_dir="benchmarks/state_space_surplus_scaling"
build_dir="build/benchmarks"
quadra_bin="$build_dir/state_space_surplus_scaling_quadra"
quadra_log="$artifact_dir/quadra_raw.txt"
tmb_log="$artifact_dir/tmb_raw.txt"
csv_path="$artifact_dir/results.csv"
plot_path="$artifact_dir/scaling_plot.png"

mkdir -p "$build_dir" "$artifact_dir"

"${CXX:-c++}" ${CXXFLAGS:--std=c++17 -O3 -DNDEBUG} \
  -I. -Iexternal/eigen -Iexternal/LBFGSpp/include \
  examples/state_space_surplus_production/benchmark_latent_tridiagonal_analytic_direct_runtime_scaled.cpp \
  -o "$quadra_bin"

# Compile the TMB DLL before measuring so compiler memory is not attributed to
# the model evaluation. The benchmark script reuses this library thereafter.
Rscript -e 'library(TMB); src <- "examples/state_space_surplus_production/tmb/state_space_surplus_tmb.cpp"; dll <- TMB::dynlib("examples/state_space_surplus_production/tmb/state_space_surplus_tmb"); if (!file.exists(dll)) TMB::compile(src, flags="-O2")'

: > "$quadra_log"
: > "$tmb_log"

measure_one() {
  local log_path="$1"
  local n="$2"
  shift 2
  local time_log
  time_log="$(mktemp)"

  echo "BEGIN n=$n" >> "$log_path"
  /usr/bin/time -l "$@" >> "$log_path" 2> "$time_log"
  cat "$time_log" >> "$log_path"

  local rss_bytes
  rss_bytes="$(awk '/maximum resident set size/ {print $1; exit}' "$time_log")"
  if [[ -z "$rss_bytes" ]]; then
    # GNU time reports KiB; normalize it to bytes.
    local rss_kib
    rss_kib="$(awk -F: '/Maximum resident set size/ {gsub(/^[[:space:]]+/, "", $2); print $2; exit}' "$time_log")"
    rss_bytes="$((rss_kib * 1024))"
  fi
  echo "PEAK_RSS_BYTES n=$n value=$rss_bytes" >> "$log_path"
  echo "END n=$n" >> "$log_path"
  rm -f "$time_log"
}

IFS=',' read -r -a sizes <<< "$lengths"
for n in "${sizes[@]}"; do
  measure_one "$quadra_log" "$n" "$quadra_bin" "$reps" "$n"
  measure_one "$tmb_log" "$n" Rscript \
    examples/state_space_surplus_production/tmb/benchmark_scaled_fixed_theta_tmb.R \
    "$reps" "$n"
done

# R's human-readable header contains a trailing blank; keep generated artifacts
# diff-clean without changing any numeric output.
perl -pi -e 's/[ \t]+$//' "$quadra_log" "$tmb_log"

python3 "$artifact_dir/normalize_results.py" \
  "$quadra_log" "$tmb_log" "$csv_path"

if ! python3 "$artifact_dir/make_scaling_plot.py" "$csv_path" "$plot_path"; then
  Rscript "$artifact_dir/make_scaling_plot.R" "$csv_path" "$plot_path"
fi

echo "CSV:  $csv_path"
echo "Plot: $plot_path"
