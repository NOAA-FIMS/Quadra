#!/usr/bin/env bash
set -euo pipefail

./inspect_obs_data_layering.sh

echo
echo "== Rebuild O3 Pollock showcase after moving Obs to data layer =="
./run_pollock_driver_showcase_report.sh
