#!/usr/bin/env bash
set -euo pipefail

./inspect_core_diagnostics_api.sh

echo
echo "== Rebuilding Pollock polished report =="
./run_pollock_polished_report.sh
