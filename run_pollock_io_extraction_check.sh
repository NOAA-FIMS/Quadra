#!/usr/bin/env bash
set -euo pipefail

./inspect_pollock_io_extraction.sh

echo
echo "== Rebuild O3 Pollock showcase after data I/O extraction =="
./run_pollock_driver_showcase_report.sh
