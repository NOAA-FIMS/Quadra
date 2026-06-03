#!/usr/bin/env bash
set -euo pipefail

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/examples

set -x
"${CXX}" ${CXXFLAGS} \
  -Iexternal/Eigen \
  -Iexternal/LBFGSpp/include \
  -Iexamples/surplus_production \
  examples/surplus_production/fit_surplus_production_lbfgs.cpp \
  -o build/examples/fit_surplus_production_lbfgs

./build/examples/fit_surplus_production_lbfgs
