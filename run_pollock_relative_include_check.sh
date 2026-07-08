#!/usr/bin/env bash
set -euo pipefail

./inspect_pollock_relative_includes.sh

echo
echo "== Rebuild O3 Pollock showcase after relative include fix =="
./run_pollock_driver_showcase_report.sh
