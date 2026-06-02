#!/usr/bin/env bash
set -euo pipefail

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2}"

EIGEN_INCLUDE=""
if [[ -d external/Eigen ]]; then
  EIGEN_INCLUDE="-Iexternal/Eigen"
elif [[ -d core/eigen ]]; then
  EIGEN_INCLUDE="-Icore/eigen"
fi

mkdir -p build/tests

set -x
"${CXX}" ${CXXFLAGS} ${EIGEN_INCLUDE} -I. \
  tests/test_exact_gradient_workspace.cpp \
  -o build/tests/test_exact_gradient_workspace

./build/tests/test_exact_gradient_workspace
