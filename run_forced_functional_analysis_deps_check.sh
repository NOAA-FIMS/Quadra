#!/usr/bin/env bash
set -euo pipefail

./inspect_forced_functional_analysis_deps.sh

echo
echo "== Rebuild O3 Pollock showcase after forced dependency move =="
./run_pollock_driver_showcase_report.sh
