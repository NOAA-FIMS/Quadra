#!/usr/bin/env bash
set -euo pipefail

mkdir -p examples/NMFS/sefsc_red_snapper/outputs

c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  -Icore \
  -o examples/NMFS/sefsc_red_snapper/quadra/red_snapper_age_structured \
  examples/NMFS/sefsc_red_snapper/quadra/red_snapper_age_structured.cpp

./examples/NMFS/sefsc_red_snapper/quadra/red_snapper_age_structured
