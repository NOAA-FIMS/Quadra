#!/usr/bin/env bash
set -euo pipefail

LENGTHS="${1:-500,1000}"
AGES="${2:-10}"

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/examples

set -x
"${CXX}" ${CXXFLAGS} \
  -Iexternal/Eigen \
  examples/age_structured_recruitment/benchmark_age_structured_no_plus_flat_band_alpha_trace.cpp \
  -o build/examples/benchmark_age_structured_no_plus_flat_band_alpha_trace

./build/examples/benchmark_age_structured_no_plus_flat_band_alpha_trace 1 "$LENGTHS" "$AGES"
