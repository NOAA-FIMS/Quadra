#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/tests

set -x
c++ -std=c++17 -O2 -DNDEBUG -g \
  -Iexternal/Eigen \
  -I. \
  tests/test_native_structured_value_update.cpp \
  -o build/tests/test_native_structured_value_update

./build/tests/test_native_structured_value_update
