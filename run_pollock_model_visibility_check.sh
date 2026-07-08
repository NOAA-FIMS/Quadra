#!/usr/bin/env bash
set -euo pipefail

./inspect_pollock_model_visibility.sh

echo
echo "== Rebuild O3 Pollock showcase after model visibility fix =="
./run_pollock_driver_showcase_report.sh
