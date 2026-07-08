#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/tests

mkdir -p build/tests

c++ -std=c++17 -O3 \
  -I. -Iexternal/eigen \
  tests/laplace/test_laplace_gradient_regression.cpp \
  examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/quadra/bigeye_adgraph_global.cpp \
  -o build/tests/test_laplace_gradient_regression

./build/tests/test_laplace_gradient_regression
