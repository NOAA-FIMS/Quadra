#!/usr/bin/env bash
set -euo pipefail

./inspect_tightened_diagnostics_api.sh

echo
echo "== Rebuild O3 clean-driver Pollock showcase with tightened diagnostics API =="
./run_pollock_driver_showcase_report.sh
