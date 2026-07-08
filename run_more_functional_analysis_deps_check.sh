#!/usr/bin/env bash
set -euo pipefail

./inspect_more_functional_analysis_deps.sh

echo
echo "== Rebuild O3 Pollock showcase after adding more FA deps =="
./run_pollock_driver_showcase_report.sh
