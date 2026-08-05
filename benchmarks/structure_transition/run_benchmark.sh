#!/usr/bin/env bash
set -euo pipefail

reps="${1:-50}"
sizes="${2:-250,1000,5000}"
out_dir="benchmarks/structure_transition"
binary="build/benchmarks/quadra_structure_transition"
csv="$out_dir/results.csv"

mkdir -p "$(dirname "$binary")" "$out_dir"
"${CXX:-c++}" ${CXXFLAGS:--std=c++17 -O3 -DNDEBUG} \
  -I. -Iexternal/eigen benchmarks/quadra_structure_transition.cpp -o "$binary"

echo "case,n,nnz,bandwidth,structure,backend,symbolic_analyses,numeric_factorizations,ms,logdet_diff,peak_rss_mib,repeated_peak_rss_mib" > "$csv"
IFS=',' read -r -a n_values <<< "$sizes"
for n in "${n_values[@]}"; do
  for kind in tridiagonal banded2 banded5 banded10 banded32 irregular; do
    repeat_log="$(mktemp)"
    time_log="$(mktemp)"
    row="$(/usr/bin/time -l "$binary" "$kind" "$n" "$reps" 2> "$repeat_log")"
    repeated_rss_bytes="$(awk '/maximum resident set size/ {print $1; exit}' "$repeat_log")"
    /usr/bin/time -l "$binary" "$kind" "$n" 1 > /dev/null 2> "$time_log"
    rss_bytes="$(awk '/maximum resident set size/ {print $1; exit}' "$time_log")"
    rm -f "$repeat_log" "$time_log"
    rss_mib="$(awk -v bytes="$rss_bytes" 'BEGIN {printf "%.6f", bytes / 1048576.0}')"
    repeated_rss_mib="$(awk -v bytes="$repeated_rss_bytes" 'BEGIN {printf "%.6f", bytes / 1048576.0}')"
    echo "$row,$rss_mib,$repeated_rss_mib" >> "$csv"
  done
done

Rscript "$out_dir/make_plot.R" "$csv" "$out_dir/structure_transition.png"
echo "CSV:  $csv"
echo "Plot: $out_dir/structure_transition.png"
