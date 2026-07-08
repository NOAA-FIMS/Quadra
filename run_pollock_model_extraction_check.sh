#!/usr/bin/env bash
set -euo pipefail

./inspect_pollock_model_extraction.sh

echo
echo "== Rebuild O3 Pollock showcase after model extraction =="
./run_pollock_driver_showcase_report.sh
