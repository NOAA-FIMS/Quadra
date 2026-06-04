#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-100}"
LENGTHS="${2:-100,500,1000,5000}"

mkdir -p build/tests benchmarks

set -x
c++ -std=c++17 -O2 -DNDEBUG -g \
  -Iexternal/Eigen \
  -I. \
  tests/benchmark_direct_structured_value_runtime_v2.cpp \
  -o build/tests/benchmark_direct_structured_value_runtime_v2

./build/tests/benchmark_direct_structured_value_runtime_v2 "$REPS" "$LENGTHS" \
  | tee "benchmarks/direct_structured_value_runtime_v2_reps${REPS}_n${LENGTHS//,/}.txt"
