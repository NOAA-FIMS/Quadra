#!/usr/bin/env bash
set -euo pipefail

./generate_bigeye_v2_caa_execution_plan.sh
python3 tools/caa/generate_assessment_cycle_from_plan.py
