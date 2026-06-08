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
  -Iexternal/LBFGSpp/include \
  examples/age_structured_recruitment/benchmark_age_structured_no_plus_analytic_banded.cpp \
  -o build/examples/benchmark_age_structured_no_plus_analytic_banded

./build/examples/benchmark_age_structured_no_plus_analytic_banded "$REPS" "$LENGTHS" "$AGES"
