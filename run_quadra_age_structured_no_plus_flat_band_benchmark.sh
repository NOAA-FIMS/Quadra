#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-10}"
LENGTHS="${2:-25,50,100,250,500,1000}"
AGES="${3:-10}"

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/examples

set -x
"${CXX}" ${CXXFLAGS} \
  -Iexternal/Eigen \
  examples/age_structured_recruitment/benchmark_age_structured_no_plus_flat_band.cpp \
  -o build/examples/benchmark_age_structured_no_plus_flat_band

./build/examples/benchmark_age_structured_no_plus_flat_band "$REPS" "$LENGTHS" "$AGES"
