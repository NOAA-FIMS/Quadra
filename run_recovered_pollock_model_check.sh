#!/usr/bin/env bash
set -euo pipefail

./inspect_recovered_pollock_model.sh

echo
echo "== Rebuild O3 Pollock showcase after recovering model header =="
./run_pollock_driver_showcase_report.sh
