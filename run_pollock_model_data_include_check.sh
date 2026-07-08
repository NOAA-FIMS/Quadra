#!/usr/bin/env bash
set -euo pipefail

./inspect_pollock_model_data_include.sh

echo
echo "== Rebuild O3 Pollock showcase after data include-depth fix =="
./run_pollock_driver_showcase_report.sh
