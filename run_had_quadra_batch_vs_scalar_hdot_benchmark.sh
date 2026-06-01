#!/usr/bin/env bash
set -euo pipefail

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O3 -DNDEBUG}"
REPS="${1:-10}"

EIGEN_INCLUDE=""
if [[ -d external/Eigen ]]; then
  EIGEN_INCLUDE="-Iexternal/Eigen"
elif [[ -d core/eigen ]]; then
  EIGEN_INCLUDE="-Icore/eigen"
fi

mkdir -p build/benchmarks

set -x
"${CXX}" ${CXXFLAGS} ${EIGEN_INCLUDE} -I. \
  benchmarks/benchmark_had_quadra_batch_vs_scalar_hdot.cpp \
  -o build/benchmarks/benchmark_had_quadra_batch_vs_scalar_hdot

./build/benchmarks/benchmark_had_quadra_batch_vs_scalar_hdot "${REPS}"
