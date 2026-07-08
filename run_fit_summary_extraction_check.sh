#!/usr/bin/env bash
set -euo pipefail

./inspect_fit_summary_extraction.sh

echo
echo "== Rebuild O3 Pollock showcase after fit-summary extraction =="
./run_pollock_driver_showcase_report.sh
