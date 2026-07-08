#!/usr/bin/env bash
set -euo pipefail

echo "== Confirm logdet-gradient sign in optimizer =="
grep -n "res.grad_x\[k\] [-+]=" -A4 -B6 core/optimizer.hpp

echo
echo "== Build/run Level 21 age-based M check =="
./run_bigeye_level21_age_based_m_check.sh
