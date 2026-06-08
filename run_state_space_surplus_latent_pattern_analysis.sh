#!/usr/bin/env bash
set -euo pipefail

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/examples

set -x
"${CXX}" ${CXXFLAGS} \
  -Iexternal/Eigen \
  -Iexternal/LBFGSpp/include \
  -I. \
  examples/state_space_surplus_production/analyze_state_space_latent_pattern.cpp \
  -o build/examples/analyze_state_space_latent_pattern

./build/examples/analyze_state_space_latent_pattern
