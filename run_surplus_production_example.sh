#!/usr/bin/env bash
set -euo pipefail

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/examples

set -x
"${CXX}" ${CXXFLAGS} -Iexamples/surplus_production   examples/surplus_production/run_surplus_production.cpp   -o build/examples/run_surplus_production

./build/examples/run_surplus_production
