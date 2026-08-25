#!/usr/bin/env bash
set -euo pipefail

reps="${1:-10}"
lengths="${2:-25,50,100,250,500,1000}"
ages="${3:-10}"

mkdir -p build/examples
"${CXX:-c++}" -std=c++17 -O2 -I. -Iexternal/eigen \
  -Iexternal/LBFGSpp/include \
  examples/age_structured_recruitment/benchmark_age_structured_recruitment.cpp \
  -o build/examples/benchmark_age_structured_recruitment

echo "Quadra benchmark"
./build/examples/benchmark_age_structured_recruitment \
  "$reps" "$lengths" "$ages"

echo
echo "TMB benchmark"
Rscript \
  examples/age_structured_recruitment/tmb/benchmark_age_structured_recruitment_tmb.R \
  "$reps" "$lengths" "$ages"
