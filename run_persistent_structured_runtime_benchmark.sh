#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-20}"
LENGTHS="${2:-100,500,1000,5000}"

mkdir -p build/tests benchmarks

set -x
c++ -std=c++17 -O2 -DNDEBUG -g \
  -Iexternal/Eigen \
  -I. \
  tests/benchmark_persistent_structured_runtime.cpp \
  -o build/tests/benchmark_persistent_structured_runtime

./build/tests/benchmark_persistent_structured_runtime "$REPS" "$LENGTHS" \
  | tee "benchmarks/persistent_structured_runtime_reps${REPS}_n${LENGTHS//,/}.txt"
