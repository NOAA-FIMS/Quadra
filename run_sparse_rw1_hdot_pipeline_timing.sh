#!/usr/bin/env bash
set -euo pipefail

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O3 -DNDEBUG -g}"
REPS="${1:-1}"

EIGEN_INCLUDE=""
if [[ -d external/Eigen ]]; then
  EIGEN_INCLUDE="-Iexternal/Eigen"
elif [[ -d core/eigen ]]; then
  EIGEN_INCLUDE="-Icore/eigen"
fi

LBFGS_INCLUDE=""
if [[ -d external/LBFGSpp/include ]]; then
  LBFGS_INCLUDE="-Iexternal/LBFGSpp/include"
elif [[ -d external/LBFGSpp ]]; then
  LBFGS_INCLUDE="-Iexternal/LBFGSpp"
fi

mkdir -p build/benchmarks

set -x
"${CXX}" ${CXXFLAGS} ${EIGEN_INCLUDE} ${LBFGS_INCLUDE} -I. \
  benchmarks/benchmark_sparse_rw1_hdot_pipeline_timing.cpp \
  -o build/benchmarks/benchmark_sparse_rw1_hdot_pipeline_timing

./build/benchmarks/benchmark_sparse_rw1_hdot_pipeline_timing "${REPS}"
