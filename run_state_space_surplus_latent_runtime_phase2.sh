#!/usr/bin/env bash
set -euo pipefail

REPS="${1:-20}"
CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/examples

set -x
"${CXX}" ${CXXFLAGS} \
  -Iexternal/Eigen \
  -Iexternal/LBFGSpp/include \
  -I. \
  examples/state_space_surplus_production/laplace_state_space_surplus_latent_runtime.cpp \
  -o build/examples/laplace_state_space_surplus_latent_runtime

./build/examples/laplace_state_space_surplus_latent_runtime "$REPS"
