#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-20}"
LENGTHS="${2:-100,500,1000,5000}"

mkdir -p build/tests benchmarks

set -x
c++ -std=c++17 -O2 -DNDEBUG -g \
  -Iexternal/Eigen \
  -I. \
  tests/benchmark_structured_value_factory.cpp \
  -o build/tests/benchmark_structured_value_factory

./build/tests/benchmark_structured_value_factory "$REPS" "$LENGTHS" \
  | tee "benchmarks/structured_value_factory_reps${REPS}_n${LENGTHS//,/}.txt"
