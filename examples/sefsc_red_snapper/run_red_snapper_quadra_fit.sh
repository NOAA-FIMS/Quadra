#!/usr/bin/env bash
set -euo pipefail

mkdir -p examples/sefsc_red_snapper/outputs

c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  -Icore \
  -Iexternal/LBFGSpp/include \
  \
  -o examples/sefsc_red_snapper/quadra/red_snapper_quadra_fit \
  examples/sefsc_red_snapper/quadra/red_snapper_quadra_fit.cpp \
  examples/sefsc_red_snapper/quadra/red_snapper_adgraph_global.cpp

./examples/sefsc_red_snapper/quadra/red_snapper_quadra_fit
