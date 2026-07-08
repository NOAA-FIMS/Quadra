#!/usr/bin/env bash
set -euo pipefail

./inspect_pollock_parameters_extraction.sh

echo
echo "== Rebuild O3 Pollock showcase after parameter extraction =="
./run_pollock_driver_showcase_report.sh
