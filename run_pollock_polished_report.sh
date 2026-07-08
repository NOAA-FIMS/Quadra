#!/usr/bin/env bash
set -euo pipefail

./run_pollock_showcase_report.sh

echo
echo "===== POLISHED POLLOCK REPORT ====="
cat examples/NMFS/afsc_walleye_pollock/outputs/walleye_pollock_analysis.md
