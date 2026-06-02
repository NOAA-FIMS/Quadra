#!/usr/bin/env bash
set -euo pipefail

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/examples

set -x
"${CXX}" ${CXXFLAGS}   examples/state_space_surplus_production/run_state_space_surplus_production.cpp   -o build/examples/run_state_space_surplus_production

./build/examples/run_state_space_surplus_production
