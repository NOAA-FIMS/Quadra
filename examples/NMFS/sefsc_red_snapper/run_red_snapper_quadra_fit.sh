#!/usr/bin/env bash
set -euo pipefail

mkdir -p examples/NMFS/sefsc_red_snapper/outputs

c++ -DQUADRA_DEBUG_FD_FINAL_GRADIENT -DQUADRA_DEBUG_LAPLACE_GRADIENT_PARTS -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  -Icore \
  -Iexternal/LBFGSpp/include \
  \
  -o examples/NMFS/sefsc_red_snapper/quadra/red_snapper_quadra_fit \
  examples/NMFS/sefsc_red_snapper/quadra/red_snapper_quadra_fit.cpp \
  examples/NMFS/sefsc_red_snapper/quadra/red_snapper_adgraph_global.cpp

./examples/NMFS/sefsc_red_snapper/quadra/red_snapper_quadra_fit
