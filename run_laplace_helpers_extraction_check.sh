#!/usr/bin/env bash
set -euo pipefail

./inspect_laplace_helpers_extraction.sh

echo
echo "== Rebuild O3 Pollock showcase after Laplace helper extraction =="
./run_pollock_driver_showcase_report.sh
