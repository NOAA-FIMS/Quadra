#!/usr/bin/env bash
set -euo pipefail

mkdir -p examples/NMFS/sefsc_red_snapper/outputs

c++ -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  -Icore \
  -o examples/NMFS/sefsc_red_snapper/quadra/evaluate_red_snapper_objective \
  examples/NMFS/sefsc_red_snapper/quadra/evaluate_red_snapper_objective.cpp

./examples/NMFS/sefsc_red_snapper/quadra/evaluate_red_snapper_objective
