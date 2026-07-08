#!/usr/bin/env bash
set -euo pipefail

./inspect_pollock_layout_normalized.sh

echo
echo "== Rebuild O3 Pollock showcase after layout cleanup =="
./run_pollock_driver_showcase_report.sh
