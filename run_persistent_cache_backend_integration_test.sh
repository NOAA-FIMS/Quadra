#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/tests

c++ -std=c++17 -O2 -DNDEBUG -g \
  -Iexternal/Eigen \
  -I. \
  tests/test_persistent_cache_backend_integration.cpp \
  -o build/tests/test_persistent_cache_backend_integration

./build/tests/test_persistent_cache_backend_integration
