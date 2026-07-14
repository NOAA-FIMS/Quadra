#!/usr/bin/env bash
set -euo pipefail
./generate_bigeye_v2_caa_ir.sh
./generate_bigeye_v2_caa_execution_plan.sh
python3 tools/caa/generate_operation_graph.py
