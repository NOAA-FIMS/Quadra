#!/usr/bin/env bash
set -euo pipefail

./inspect_pollock_constants_extraction.sh

echo
echo "== Rebuild O3 Pollock showcase after constants extraction =="
./run_pollock_driver_showcase_report.sh
