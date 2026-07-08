#!/usr/bin/env bash
set -euo pipefail

./inspect_functional_analysis_dependency_repair.sh

echo
echo "== Rebuild O3 Pollock showcase after functional-analysis dependency repair =="
./run_pollock_driver_showcase_report.sh
