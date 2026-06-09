#!/usr/bin/env bash
set -euo pipefail

echo "== Inspect existing Quadra H_uu extraction patterns =="

echo
echo "== core/autodiff.hpp ADScope section =="
sed -n '100,145p' core/autodiff.hpp || true

echo
echo "== core/laplace.hpp extract_sparse_hessian signature =="
sed -n '280,330p' core/laplace.hpp || true

echo
echo "== core/laplace.hpp first in-repo extraction usage around 480-525 =="
sed -n '470,530p' core/laplace.hpp || true

echo
echo "== core/laplace.hpp extraction usage around 545-590 =="
sed -n '545,590p' core/laplace.hpp || true

echo
echo "== core/laplace.hpp extraction usage around 700-735 =="
sed -n '700,735p' core/laplace.hpp || true

echo
echo "== core/laplace.hpp extraction usage around 780-810 =="
sed -n '780,810p' core/laplace.hpp || true

echo
echo "== constructors / independent variable helpers =="
grep -RIn "Make.*AD\\|Independent\\|ADScope scope\\|std::vector<AD> p_full\\|p_full.push_back" core examples tests | head -120 || true
