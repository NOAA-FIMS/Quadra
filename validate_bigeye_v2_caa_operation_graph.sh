#!/usr/bin/env bash
set -euo pipefail

./generate_bigeye_v2_caa_operation_graph.sh
./generate_bigeye_v2_caa_execution_plan.sh
python3 tools/caa/validate_operation_graph.py
