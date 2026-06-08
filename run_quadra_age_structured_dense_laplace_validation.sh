#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-3}"
LENGTHS="${2:-25,50}"
AGES="${3:-10}"

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/examples

set -x
"${CXX}" ${CXXFLAGS} \
  -Iexternal/Eigen \
  -Iexternal/LBFGSpp/include \
  examples/age_structured_recruitment/benchmark_age_structured_dense_laplace_validation.cpp \
  -o build/examples/benchmark_age_structured_dense_laplace_validation

./build/examples/benchmark_age_structured_dense_laplace_validation "$REPS" "$LENGTHS" "$AGES"
