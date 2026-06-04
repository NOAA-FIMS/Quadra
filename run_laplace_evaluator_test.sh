#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/tests

set -x
c++ -std=c++17 -O2 -DNDEBUG -g \
  -Iexternal/Eigen \
  -I. \
  tests/test_laplace_evaluator.cpp \
  -o build/tests/test_laplace_evaluator

./build/tests/test_laplace_evaluator
