#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$repo_root"

mkdir -p build/examples examples/NMFS/sefsc_red_snapper/outputs

"${CXX:-c++}" -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  -Icore \
  -Iexternal/LBFGSpp/include \
  -o build/examples/sefsc_red_snapper_quadra_fit \
  examples/NMFS/sefsc_red_snapper/quadra/red_snapper_quadra_fit.cpp \
  examples/NMFS/sefsc_red_snapper/quadra/red_snapper_adgraph_global.cpp

./build/examples/sefsc_red_snapper_quadra_fit
