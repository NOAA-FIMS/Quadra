#!/usr/bin/env bash
set -euo pipefail

./inspect_functional_analysis_diagnostics_extraction.sh

echo
echo "== Rebuild O3 Pollock showcase after functional-analysis diagnostics extraction =="
./run_pollock_driver_showcase_report.sh
