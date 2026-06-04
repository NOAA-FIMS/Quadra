#!/usr/bin/env bash
set -euo pipefail

mkdir -p build/tests

set -x
c++ -std=c++17 -O2 -DNDEBUG -g \
  -Iexternal/Eigen \
  -I. \
  tests/test_persistent_latent_state_runtime.cpp \
  -o build/tests/test_persistent_latent_state_runtime

./build/tests/test_persistent_latent_state_runtime
