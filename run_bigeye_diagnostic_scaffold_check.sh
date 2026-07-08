#!/usr/bin/env bash
set -euo pipefail

./inspect_bigeye_diagnostic_scaffold.sh

echo
echo "== Git status for Bigeye scaffold =="
git status --short examples/NMFS/pifsc_bigeye_tuna
