#!/usr/bin/env bash
set -euo pipefail

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/tests

set -x
"${CXX}" ${CXXFLAGS} \
  -Iexternal/Eigen \
  -I. \
  tests/test_laplace_backend_factory.cpp \
  -o build/tests/test_laplace_backend_factory

./build/tests/test_laplace_backend_factory
