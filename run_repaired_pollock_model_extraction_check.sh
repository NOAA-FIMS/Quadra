#!/usr/bin/env bash
set -euo pipefail

./inspect_repaired_pollock_model_extraction.sh

echo
echo "== Rebuild O3 Pollock showcase after repaired model extraction =="
./run_pollock_driver_showcase_report.sh
