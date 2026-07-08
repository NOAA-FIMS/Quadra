#!/usr/bin/env bash
set -euo pipefail
L21="examples/NMFS/pifsc_bigeye_tuna/level21_age_based_natural_mortality_diagnostic"

echo "== Shared layout header =="
cat "$L21/diagnostics/bigeye_level21_parameter_layout.hpp"

echo
echo "== Duplicate bare constant definitions scan =="
grep -R "constexpr int kBaseFixed\|constexpr int kMParamOffset\|constexpr int kLonglineSelOffset" -n \
  "$L21/diagnostics" "$L21/reports" "$L21/objective" || true

echo
echo "== Includes layout header =="
grep -R "bigeye_level21_parameter_layout.hpp\|using namespace level21_layout" -n \
  "$L21/diagnostics" "$L21/reports" || true
