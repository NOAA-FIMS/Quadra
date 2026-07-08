#!/usr/bin/env bash
set -euo pipefail

./inspect_fixed_effect_diagnostics_extraction.sh

echo
echo "== Rebuild O3 Pollock showcase after fixed-effect diagnostics extraction =="
./run_pollock_driver_showcase_report.sh
