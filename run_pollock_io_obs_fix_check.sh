#!/usr/bin/env bash
set -euo pipefail

./inspect_pollock_io_obs_fix.sh

echo
echo "== Rebuild O3 Pollock showcase after I/O Obs include fix =="
./run_pollock_driver_showcase_report.sh
