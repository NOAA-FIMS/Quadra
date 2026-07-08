#!/usr/bin/env bash
set -euo pipefail

./inspect_bigeye_scientific_workflow_docs.sh

echo
echo "== Git status for Bigeye workflow docs =="
git status --short examples/NMFS/pifsc_bigeye_tuna/workflow
