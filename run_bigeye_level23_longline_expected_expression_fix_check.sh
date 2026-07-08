#!/usr/bin/env bash
set -euo pipefail

ROOT="examples/NMFS/pifsc_bigeye_tuna/level23_longline_selectivity_smoothness_scan"
LL="$ROOT/diagnostics/bigeye_longline_prediction_decomposition.hpp"

echo "== Longline BH context =="
grep -n "const double r0\|default_weight_at_age\|default_maturity_at_age\|const double phi0\|spawning_biomass_from_numbers\|expected_recruitment\|next\[0\]" "$LL"

echo
echo "== Build/run Level 23 BH check =="
./run_bigeye_level23_bh_diag_compile_fix_check.sh
