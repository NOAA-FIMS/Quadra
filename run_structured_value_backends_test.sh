#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/tests

set -x
c++ -std=c++17 -O2 -DNDEBUG -g \
  -Iexternal/Eigen \
  -I. \
  tests/test_structured_value_backends.cpp \
  -o build/tests/test_structured_value_backends

./build/tests/test_structured_value_backends
