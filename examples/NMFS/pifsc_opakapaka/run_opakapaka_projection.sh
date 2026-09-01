#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$repo_root"

mkdir -p build/examples examples/NMFS/pifsc_opakapaka/outputs

"${CXX:-c++}" -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  -Iexternal/LBFGSpp/include \
  examples/NMFS/pifsc_opakapaka/quadra/opakapaka.cpp \
  examples/NMFS/pifsc_opakapaka/quadra/opakapaka_adgraph_global.cpp \
  -o build/examples/pifsc_opakapaka

./build/examples/pifsc_opakapaka
