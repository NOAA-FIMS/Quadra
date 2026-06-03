#!/usr/bin/env bash
set -euo pipefail

N="${1:-25}"
AGES="${2:-10}"

CXX="${CXX:-c++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -DNDEBUG -g}"

mkdir -p build/examples

set -x
"${CXX}" ${CXXFLAGS} \
  -Iexternal/Eigen \
  -Iexternal/LBFGSpp/include \
  -Iexamples/age_structured_recruitment \
  examples/age_structured_recruitment/check_age_structured_no_plus_derivatives.cpp \
  -o build/examples/check_age_structured_no_plus_derivatives

./build/examples/check_age_structured_no_plus_derivatives "$N" "$AGES"
