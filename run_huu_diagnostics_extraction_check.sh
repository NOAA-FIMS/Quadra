#!/usr/bin/env bash
set -euo pipefail

./inspect_huu_diagnostics_extraction.sh

echo
echo "== Rebuild O3 Pollock showcase after Huu diagnostics extraction =="
./run_pollock_driver_showcase_report.sh
