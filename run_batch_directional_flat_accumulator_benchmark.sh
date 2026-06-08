#!/usr/bin/env bash
set -euo pipefail

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O3 -DNDEBUG -g}"
REPS="${1:-100}"

EIGEN_INCLUDE=""
if [[ -d external/Eigen ]]; then
  EIGEN_INCLUDE="-Iexternal/Eigen"
elif [[ -d core/eigen ]]; then
  EIGEN_INCLUDE="-Icore/eigen"
fi

mkdir -p build/benchmarks

set -x
"${CXX}" ${CXXFLAGS} ${EIGEN_INCLUDE} -I. \
  benchmarks/benchmark_batch_directional_flat_accumulator.cpp \
  -o build/benchmarks/benchmark_batch_directional_flat_accumulator

./build/benchmarks/benchmark_batch_directional_flat_accumulator "${REPS}"
