#!/usr/bin/env bash
set -euo pipefail

OUT=tests/big_laplace_convergence_contract.out
BUILD_LOG=tests/big_laplace_convergence_contract.build.log
EXE=examples/big/catch_at_age_laplace

CXX="${CXX:-g++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include}"

echo "Compiling $EXE"
if ! "$CXX" $CXXFLAGS -o "$EXE" examples/big/catch_at_age_laplace.cpp > "$BUILD_LOG" 2>&1; then
  echo "ERROR: contract compile failed"
  cat "$BUILD_LOG"
  exit 2
fi

echo "Running $EXE"
if ! "$EXE" > "$OUT" 2>&1; then
  echo "ERROR: big Laplace example failed"
  cat "$OUT"
  exit 2
fi

python3 - <<'PY'
from pathlib import Path
import re, sys, math

txt = Path("tests/big_laplace_convergence_contract.out").read_text()

penalty_count = txt.count("returning penalty")
if penalty_count:
    print(f"NOTE: Laplace penalty path occurred during trial evaluations: {penalty_count}")

fit = re.search(r"fit value:\s*([-+0-9.eE]+)", txt)
grad = re.search(r"optimizer-reported fixed gradient norm:\s*([-+0-9.eE]+)", txt)
direct = re.search(r"direct_full_minus_reported:\s*([-+0-9.eE]+)", txt)

if not (fit and grad and direct):
    print("FAIL: could not parse contract diagnostics")
    sys.exit(1)

fit = float(fit.group(1))
grad = float(grad.group(1))
direct = abs(float(direct.group(1)))

ok = (
    math.isfinite(fit)
    and -500.0 < fit < 0.0
    and grad < 1.0e-3
    and direct < 1.0e-8
)

if not ok:
    print("FAIL: big Laplace convergence contract failed")
    print(f"  fit value: {fit}")
    print(f"  fixed gradient norm: {grad}")
    print(f"  direct_full_minus_reported abs: {direct}")
    sys.exit(1)

print("PASS: big Laplace black-box convergence contract satisfied")
print(f"  fit value: {fit}")
print(f"  fixed gradient norm: {grad}")
print(f"  direct_full_minus_reported abs: {direct}")
PY
