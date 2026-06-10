#!/usr/bin/env bash
set -euo pipefail

mkdir -p examples/sefsc_red_snapper/outputs

c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  -Icore \
  -o examples/sefsc_red_snapper/quadra/red_snapper_level0 \
  examples/sefsc_red_snapper/quadra/red_snapper_level0.cpp

./examples/sefsc_red_snapper/quadra/red_snapper_level0
