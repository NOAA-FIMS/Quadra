#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/examples examples/opakapaka_projection/outputs

set -x
c++ -std=c++17 -O2 -DNDEBUG -g \
  -Iexternal/Eigen \
  -Iexternal/LBFGSpp/include \
  -I. \
  examples/opakapaka_projection/opakapaka_synthetic_fit_projection.cpp \
  -o build/examples/opakapaka_synthetic_fit_projection

./build/examples/opakapaka_synthetic_fit_projection
