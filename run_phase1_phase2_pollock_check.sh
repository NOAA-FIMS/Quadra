#!/usr/bin/env bash
set -euo pipefail

./inspect_phase1_phase2_pollock.sh

echo
echo "== Building/running clean Pollock driver =="
./run_pollock_driver_showcase_report.sh
