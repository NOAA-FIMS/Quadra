#!/usr/bin/env bash
set -euo pipefail

./inspect_final_driver_polish.sh

echo
echo "== Rebuild O3 Pollock showcase after final driver polish =="
./run_pollock_driver_showcase_report.sh
