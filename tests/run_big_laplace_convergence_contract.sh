#!/usr/bin/env bash
set -euo pipefail
set -x

OUT=tests/big_laplace_convergence_contract.out
BUILD_LOG=tests/big_laplace_convergence_contract.build.log

if ! make run-big-examples CXXFLAGS="-std=c++17 -O3 -flto -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include" > "$BUILD_LOG" 2>&1; then
  echo "ERROR: make run-big-examples failed"
  cat "$BUILD_LOG"
  exit 2
fi

cat "$BUILD_LOG" > "$OUT"

grep -q "PASS: big Laplace black-box convergence contract satisfied" "$OUT" || {
  echo "ERROR: PASS marker missing"
  cat "$OUT"
  exit 2
}
