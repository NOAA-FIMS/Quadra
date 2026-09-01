#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$repo_root"

mkdir -p build/examples examples/NMFS/sefsc_red_snapper/outputs

"${CXX:-c++}" -std=c++17 -O3 \
  -I. \
  -Iexternal/eigen \
  -Icore \
  -o build/examples/sefsc_red_snapper_objective \
  examples/NMFS/sefsc_red_snapper/quadra/evaluate_red_snapper_objective.cpp

./build/examples/sefsc_red_snapper_objective
