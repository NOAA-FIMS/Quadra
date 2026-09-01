#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
cd "$repo_root"

mkdir -p build/examples
"${CXX:-c++}" -std=c++17 -O3 -I. -Iexternal/eigen \
  examples/NMFS/afsc_walleye_pollock/quadra/walleye_pollock.cpp \
  examples/NMFS/afsc_walleye_pollock/quadra/walleye_pollock_adgraph_global.cpp \
  -o build/examples/afsc_walleye_pollock
./build/examples/afsc_walleye_pollock
