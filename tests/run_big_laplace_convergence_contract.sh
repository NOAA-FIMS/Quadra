#!/usr/bin/env bash
set -x
set -euo pipefail

OUT="tests/big_laplace_convergence_contract.out"

make run-big-examples \
  CXXFLAGS='-std=c++17 -O3 -flto -I. -I./external/eigen -I./external/had -I./external/LBFGSpp/include' \
  > "${OUT}" 2>&1

python3 - <<'PY'
from pathlib import Path
import math
import re
import sys

out_path = Path("tests/big_laplace_convergence_contract.out")
text = out_path.read_text()

def fail(msg):
    print(f"FAIL: {msg}")
    print(f"See: {out_path}")
    sys.exit(1)

def grab(pattern, name):
    m = re.search(pattern, text)
    if not m:
        fail(f"missing {name}")
    try:
        return float(m.group(1))
    except ValueError:
        fail(f"could not parse {name}: {m.group(1)!r}")

if "Full direct model objective check" not in text:
    fail("missing full direct model objective check block")

if "L-BFGS minimize status report" not in text:
    fail("missing L-BFGS minimize status report block")

direct_minus = grab(r"direct_full_minus_reported:\s*([-+0-9.eE]+)", "direct_full_minus_reported")
random_grad = grab(r"random gradient norm:\s*([-+0-9.eE]+)", "random gradient norm")
fit_grad = grab(r"optimizer-reported fixed gradient norm:\s*([-+0-9.eE]+)", "optimizer fixed gradient norm")
fit_value = grab(r"fit value:\s*([-+0-9.eE]+)", "fit value")
laplace_objective = grab(r"Laplace objective:\s*([-+0-9.eE]+)", "Laplace objective")

if not math.isfinite(fit_value):
    fail("fit value is not finite")
if not math.isfinite(laplace_objective):
    fail("Laplace objective is not finite")
if abs(direct_minus) > 1e-8:
    fail(f"direct_full_minus_reported too large: {direct_minus}")
if random_grad > 1e-6:
    fail(f"random-effect mode gradient too large: {random_grad}")

# This is intentionally relaxed. We learned the example is scientifically stable
# well before machine-tight outer gradients, and LBFGSpp iteration accounting can
# be misleading when line searches are involved.
if fit_grad > 1e-4:
    fail(f"fixed-effect gradient exceeds relaxed contract: {fit_grad}")

print("PASS: big Laplace black-box convergence contract satisfied")
print(f"  fit value: {fit_value}")
print(f"  Laplace objective: {laplace_objective}")
print(f"  fixed gradient norm: {fit_grad}")
print(f"  random gradient norm: {random_grad}")
print(f"  direct_full_minus_reported: {direct_minus}")
PY
