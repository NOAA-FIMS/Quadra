#!/usr/bin/env bash
set -euo pipefail

./generate_bigeye_v2_caa_ir.sh
python3 tools/caa/generate_assessment_cycle_from_ir.py
