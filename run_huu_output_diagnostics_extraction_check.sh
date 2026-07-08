#!/usr/bin/env bash
set -euo pipefail

./inspect_huu_output_diagnostics_extraction.sh

echo
echo "== Rebuild O3 Pollock showcase after Huu output extraction =="
./run_pollock_driver_showcase_report.sh
