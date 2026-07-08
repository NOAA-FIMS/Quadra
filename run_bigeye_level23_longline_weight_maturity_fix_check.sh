#!/usr/bin/env bash
set -euo pipefail

ROOT="examples/NMFS/pifsc_bigeye_tuna/level23_longline_selectivity_smoothness_scan"
LL="$ROOT/diagnostics/bigeye_longline_prediction_decomposition.hpp"

echo "== Longline diagnostic BH local context =="
grep -n "default_weight_at_age\|default_maturity_at_age\|spawning_biomass_from_numbers\|const double phi0\|expected_recruitment\|next\[0\]" "$LL"

echo
echo "== Build/run Level 23 BH diagnostic sync check =="
./run_bigeye_level23_bh_diag_compile_fix_check.sh
