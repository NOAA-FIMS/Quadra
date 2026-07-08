#!/usr/bin/env bash
set -euo pipefail

echo "== Gradient parts FD audit insertion =="
grep -n "level21_gradient_parts_fd_audit_v1\|Gradient Parts FD Audit\|total_cosine\|joint_cosine\|logdet_cosine" \
  -A12 -B8 \
  examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic/quadra/bigeye_level21_age_based_natural_mortality_diagnostic.cpp

echo
echo "== Run Level 21 with gradient parts FD audit =="
./run_bigeye_level21_age_based_m_check.sh
